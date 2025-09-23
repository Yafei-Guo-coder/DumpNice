#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
脚本4：方案一分析 - 使用scGLUE推断生成的feather文件修剪网络
功能：实现方案一的完整GRN构建流程，包括GRN草稿生成和使用scGLUE推断的feather文件进行网络修剪
特点：始终使用HVG（高可变基因）和TF（转录因子）的并集作为基因集
参数说明：
  --top_genes：控制使用多少高可变基因（默认0，表示使用所有高可变基因）
  --rank_threshold：数值越小表示筛选越严格，限制更多的TF-DNA结合
使用方法：python 04_grn_analysis_scheme1.py --rna_h5ad <rna_h5ad_file> --tf_file <tf_file> --integration_dir <integration_results_dir> --output_dir <output_directory>
"""

import os
import argparse
import subprocess
import scanpy as sc
import pandas as pd
import numpy as np
import loompy as lp
import datetime
import logging
import platform
import warnings

# 忽略pkg_resources弃用警告
warnings.filterwarnings("ignore", category=DeprecationWarning, module="pkg_resources")

# 设置日志
def setup_logging(log_file, log_level=logging.INFO):
    """设置日志配置"""
    # 创建日志器
    logger = logging.getLogger()
    logger.setLevel(log_level)
    
    # 避免重复添加处理器
    if logger.handlers:
        logger.handlers.clear()
    
    # 创建文件处理器 - 记录所有级别的日志
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    
    # 创建控制台处理器 - 记录INFO及以上级别的日志
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    
    # 设置详细的日志格式
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    console_handler.setFormatter(console_formatter)
    
    # 添加处理器到日志器
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

# 解析命令行参数
def parse_arguments():
    parser = argparse.ArgumentParser(description='GRN Analysis Scheme 1: scGLUE-based GRN Construction and Analysis with Efficiency')
    
    # 输入文件参数
    parser.add_argument('--rna_h5ad', type=str, required=True, help='Path to RNA h5ad file')
    parser.add_argument('--tf_file', type=str, required=True, help='Path to transcription factors file')
    parser.add_argument('--integration_dir', type=str, required=True, help='Directory containing integration results from scGLUE')
    
    # 输出参数
    parser.add_argument('--output_dir', type=str, default='./results/scheme1', help='Output directory for results')
    parser.add_argument('--output_prefix', type=str, default='scheme1_', help='Prefix for output files')
    
    # 分析参数
    parser.add_argument('--num_workers', type=int, default=4, help='Number of workers for parallel processing')
    parser.add_argument('--rank_threshold', type=int, default=2000, help='Rank threshold for regulatory regions (数值越小表示筛选越严格，限制更多的TF-DNA结合)')
    parser.add_argument('--min_genes', type=int, default=20, help='Minimum number of genes per regulon')
    parser.add_argument('--top_genes', type=int, default=0, help='Number of top highly variable genes to use (0 to use all genes)，默认使用所有HVG')
    
    # 可选参数
    parser.add_argument('--feather_dir', type=str, default=None, help='Alternative directory for feather files (if different from integration_dir)')
    
    return parser.parse_args()

# 运行命令的辅助函数
def run_command(cmd, logger=None):
    """在当前环境中运行命令并记录详细信息"""
    # 构建完整的命令字符串用于记录
    cmd_str = ' '.join(cmd)
    
    if logger:
        logger.info(f"执行命令: {cmd_str}")
        logger.debug(f"命令参数详情: {cmd}")
    else:
        print(f"Running command: {cmd_str}")
    
    try:
        # 启动进程并实时获取输出
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1  # 行缓冲
        )
        
        # 实时读取并记录标准输出
        stdout_lines = []
        for line in process.stdout:
            line = line.rstrip()
            stdout_lines.append(line)
            if logger:
                # 根据输出内容判断日志级别
                if any(keyword in line.lower() for keyword in ['error', 'exception', 'failed']):
                    logger.error(line)
                elif 'warning' in line.lower():
                    logger.warning(line)
                elif 'info' in line.lower() or 'success' in line.lower():
                    logger.info(line)
                else:
                    logger.debug(line)
            else:
                print(f"Command output: {line}")
        
        # 读取并记录标准错误
        stderr = process.stderr.read().strip()
        if stderr:
            if logger:
                logger.error(f"命令错误输出: {stderr}")
            else:
                print(f"Warning: {stderr}")
        
        # 等待进程完成并获取返回码
        process.wait()
        
        if process.returncode != 0:
            error_msg = f"命令执行失败，返回码: {process.returncode}"
            if logger:
                logger.error(error_msg)
                logger.error(f"标准输出完整内容: {' '.join(stdout_lines)}")
                logger.error(f"标准错误完整内容: {stderr}")
            else:
                print(error_msg)
            raise subprocess.CalledProcessError(process.returncode, cmd)
        
        if logger:
            logger.debug(f"命令执行成功，返回码: {process.returncode}")
            logger.debug(f"标准输出行数: {len(stdout_lines)}")
            if stdout_lines:
                logger.debug(f"标准输出前5行: {' '.join(stdout_lines[:5])}")
        
        return 0
    except Exception as e:
        if logger:
            logger.exception(f"运行命令时发生异常: {e}")
            logger.error(f"异常类型: {type(e).__name__}")
        else:
            print(f"Command failed with error: {e}")
        raise

# 生成GRN草稿
def generate_grn_draft(rna_h5ad, tf_file, output_dir, output_prefix, num_workers, top_genes=None, logger=None):
    """生成GRN草稿网络"""
    # 如果没有提供logger，使用默认的print函数
    if logger is None:
        log_info = print
    else:
        log_info = logger.info
    
    log_info("===== Step 1: Generate GRN Draft ======")
    
    # 加载snRNA数据
    log_info(f"Loading snRNA data from {rna_h5ad}...")
    if logger: logger.debug(f"Data file size: {os.path.getsize(rna_h5ad)/1024/1024:.2f} MB")
    rna = sc.read_h5ad(rna_h5ad)
    if logger: logger.debug(f"Loaded RNA data with {rna.n_obs} cells and {rna.n_vars} genes")
    
    # 加载TF列表
    log_info(f"Loading TF list from {tf_file}...")
    tfs = pd.read_csv(tf_file, header=None, names=["gene_id"])["gene_id"].values
    if logger: logger.debug(f"Loaded {len(tfs)} TFs from file")
    
    # 筛选表达的TF
    log_info("Filtering expressed TFs...")
    expressed_tfs = pd.Index(tfs).intersection(rna.var_names)
    log_info(f"Number of expressed TFs: {len(expressed_tfs)}")
    if logger: logger.debug(f"TFs found in data: {len(expressed_tfs)}/{len(tfs)} ({len(expressed_tfs)/len(tfs)*100:.1f}%)")
    
    # 保存TF列表
    tfs_output = os.path.join(output_dir, f"{output_prefix}tfs.txt")
    np.savetxt(tfs_output, expressed_tfs, fmt="%s")
    log_info(f"Expressed TFs saved to {tfs_output}")
    
    # 计算高可变基因
    log_info(f"Calculating highly variable genes...")
    
    # 预处理数据：处理可能存在的无穷大值和NaN值
    log_info(f"Preprocessing data to handle infinity values and NaNs...")
    
    # 创建数据副本以避免修改原始数据
    rna_processed = rna.copy()
    
    # 检查数据是否包含无穷大值或NaN值
    has_inf = np.isinf(rna_processed.X.data).any()
    has_nan = np.isnan(rna_processed.X.data).any()
    
    if has_inf or has_nan:
        log_info(f"Found {has_inf and 'infinity' or ''}{has_inf and has_nan and ' and ' or ''}{has_nan and 'NaN' or ''} values in data, will preprocess...")
        
        # 将无穷大值和NaN值替换为0
        if isinstance(rna_processed.X, np.ndarray):
            rna_processed.X[np.isinf(rna_processed.X)] = 0
            rna_processed.X[np.isnan(rna_processed.X)] = 0
        else:
            # 对于稀疏矩阵
            data = rna_processed.X.data
            data[np.isinf(data)] = 0
            data[np.isnan(data)] = 0
    
    if top_genes and top_genes > 0 and top_genes < rna.n_vars:
        # 使用指定数量的高可变基因
        log_info(f"Using top {top_genes} highly variable genes")
        # 使用cell_ranger flavor可以更好地处理异常值
        sc.pp.highly_variable_genes(rna_processed, n_top_genes=top_genes, flavor="cell_ranger")
    else:
        # 默认使用所有高可变基因
        log_info(f"Using all highly variable genes (default)")
        sc.pp.highly_variable_genes(rna_processed, flavor="cell_ranger", n_top_genes=None)
    
    # 创建HVG掩码
    hvg_mask = rna_processed.var['highly_variable']
    
    # 确保所有TFs都被保留在基因集中（HVG和TF的并集）
    log_info(f"Ensuring all expressed TFs are included in the gene set...")
    for tf in expressed_tfs:
        if tf in rna_processed.var_names:
            hvg_mask.loc[tf] = True
    
    # 应用筛选，使用HVG和TF的并集
    rna_subset = rna[:, hvg_mask].copy()
    log_info(f"Selected {rna_subset.n_vars} genes (including {len(expressed_tfs)} TFs) for GRN construction")
    if logger: 
        logger.debug(f"RNA subset shape after gene selection: {rna_subset.shape}")
        logger.debug(f"Number of TFs in subset: {sum(1 for tf in expressed_tfs if tf in rna_subset.var_names)}")
        logger.debug(f"Number of highly variable genes: {sum(rna_processed.var['highly_variable'])}")
    
    # 创建loom文件
    log_info("Creating loom file for pySCENIC...")
    # 使用安全的方式修改obs，避免ImplicitModificationWarning
    cells_values = rna_subset.obs_names
    # 创建新的DataFrame来避免直接修改
    obs_df = rna_subset.obs.copy()
    obs_df['cells'] = cells_values
    rna_subset.obs = obs_df
    
    loom_output = os.path.join(output_dir, f"{output_prefix}rna.loom")
    
    # 使用loompy创建loom文件
    row_attrs = {
        "Gene": np.array(rna_subset.var_names),
        "Accession": np.array(rna_subset.var_names)
    }
    col_attrs = {
        "CellID": np.array(rna_subset.obs_names),
        "cells": np.array(rna_subset.obs_names)
    }
    
    # 创建loom文件
    lp.create(loom_output, rna_subset.X.transpose(), row_attrs, col_attrs)
    log_info(f"Loom file created with {rna_subset.shape[1]} genes and {rna_subset.shape[0]} cells")
    log_info(f"Loom file saved to {loom_output}")
    if logger: logger.debug(f"Loom file size: {os.path.getsize(loom_output)/1024/1024:.2f} MB")
    
    # 执行GRN构建（pyscenic grn）
    log_info("Running pyscenic grn to build gene regulatory network draft...")
    grn_output = os.path.join(output_dir, f"{output_prefix}draft_grn.csv")
    
    # 构建pyscenic grn命令
    grn_cmd = [
        "pyscenic", "grn",
        loom_output, tfs_output,
        "-o", grn_output,
        "--seed", "0",
        "--num_workers", str(num_workers),
        "--cell_id_attribute", "cells",
        "--gene_attribute", "Gene"
    ]
    
    # 运行命令
    if logger: logger.debug(f"Running command: {' '.join(grn_cmd)}")
    try:
        # 直接在当前环境中运行pyscenic
        subprocess.run(grn_cmd, check=True)
    except Exception as e:
        if logger:
            logger.error(f"pyscenic command failed: {e}")
        else:
            print(f"pyscenic command failed: {e}")
        raise
    
    log_info(f"GRN draft construction completed, output saved to {grn_output}")
    if logger:
        if os.path.exists(grn_output):
            logger.debug(f"GRN file size: {os.path.getsize(grn_output)/1024/1024:.2f} MB")
            # 检查文件内容行数
            with open(grn_output, 'r') as f:
                line_count = sum(1 for _ in f)
            logger.debug(f"GRN file contains {line_count} lines")
    
    return loom_output, grn_output

# 使用scGLUE推断的feather文件修剪网络
def prune_network_with_glue(grn_file, loom_file, integration_dir, output_dir, output_prefix, num_workers, rank_threshold, min_genes, feather_dir=None, logger=None):
    """使用scGLUE推断的feather文件修剪网络"""
    # 如果没有提供logger，使用默认的print函数
    if logger is None:
        log_info = print
    else:
        log_info = logger.info
    
    log_info("\n===== Step 2: Prune Network with scGLUE-inferred Feather Files ======")
    
    # 决定在哪里查找feather文件
    search_dir = feather_dir if feather_dir is not None else integration_dir
    
    # 递归查找指定目录及其子目录中的feather文件
    log_info(f"Searching for feather files in {search_dir} and its subdirectories...")
    if logger: 
        logger.debug(f"Integration directory: {integration_dir}")
        if feather_dir:
            logger.debug(f"Feather directory (explicitly provided): {feather_dir}")
    
    feather_files = []
    
    # 使用os.walk递归查找所有feather文件
    for root, dirs, files in os.walk(search_dir):
        for file in files:
            if file.endswith('.feather'):
                feather_path = os.path.join(root, file)
                feather_files.append(feather_path)
                log_info(f"Found feather file: {feather_path}")
                if logger: logger.debug(f"Feather file size: {os.path.getsize(feather_path)/1024/1024:.2f} MB")
    
    if not feather_files:
        error_msg = f"No feather files found in {search_dir}"
        if logger: logger.error(error_msg)
        raise FileNotFoundError(error_msg)
    
    log_info(f"Total feather files found: {len(feather_files)}")
    if logger: logger.debug(f"Feather files list: {feather_files}")
    
    # 查找注释文件（优先在feather_dir中查找，然后在integration_dir中查找）
    annotation_file = None
    search_dirs = []
    
    # 如果提供了feather_dir，先在feather_dir中查找
    if feather_dir:
        search_dirs.append(feather_dir)
        log_info(f"Searching for annotation file in feather directory: {feather_dir}")
    
    # 始终在integration_dir中查找
    search_dirs.append(integration_dir)
    log_info(f"Searching for annotation file in integration directory: {integration_dir}")
    
    # 在指定的目录中查找注释文件
    for search_path in search_dirs:
        # 先尝试查找以.tsv结尾且包含'annotation'或'ctx'的文件
        if os.path.exists(search_path):
            for file in os.listdir(search_path):
                if file.endswith('.tsv') and ('annotation' in file.lower() or 'ctx' in file.lower()):
                    annotation_file = os.path.join(search_path, file)
                    log_info(f"Found annotation file: {annotation_file}")
                    if logger: logger.debug(f"Annotation file size: {os.path.getsize(annotation_file)/1024:.2f} KB")
                    break
        
        if annotation_file:
            break
        
        # 如果没找到，尝试使用默认命名
        default_annotation = os.path.join(search_path, "ctx_annotation.tsv")
        if os.path.exists(default_annotation):
            annotation_file = default_annotation
            log_info(f"Using default annotation file: {annotation_file}")
            break
    
    # 如果都没找到，抛出错误
    if not annotation_file:
        error_msg = f"No annotation file found in specified directories: {', '.join(search_dirs)}"
        if logger: logger.error(error_msg)
        raise FileNotFoundError(error_msg)
    
    # 执行调控子修剪（pyscenic ctx）
    log_info("Running pyscenic ctx to prune regulatory network using scGLUE-inferred data...")
    ctx_output = os.path.join(output_dir, f"{output_prefix}pruned_grn.csv")
    
    # 构建pyscenic ctx命令
    ctx_cmd = [
        "pyscenic", "ctx",
        grn_file
    ] + feather_files + [
        "--annotations_fname", annotation_file,
        "--expression_mtx_fname", loom_file,
        "--output", ctx_output,
        "--rank_threshold", str(rank_threshold),
        "--min_genes", str(min_genes),
        "--num_workers", str(num_workers),
        "--cell_id_attribute", "cells",
        "--gene_attribute", "Gene"
    ]
    
    # 运行命令
    try:
        run_command(ctx_cmd, logger)
    except subprocess.CalledProcessError as e:
        error_msg = f"pyscenic ctx command failed: {e}"
        if logger: logger.error(error_msg)
        else:
            print(error_msg)
        raise
    
    log_info(f"Regulatory network pruning completed, output saved to {ctx_output}")
    if logger: 
        if os.path.exists(ctx_output):
            logger.debug(f"Pruned GRN file size: {os.path.getsize(ctx_output)/1024/1024:.2f} MB")
    
    # 计算调控子活性（pyscenic aucell）
    log_info("Running pyscenic aucell to calculate regulon activity...")
    aucell_output = os.path.join(output_dir, f"{output_prefix}regulon_activity.csv")
    
    # 构建pyscenic aucell命令
    aucell_cmd = [
        "pyscenic", "aucell",
        loom_file,
        ctx_output,
        "--output", aucell_output,
        "--num_workers", str(num_workers)
    ]
    
    # 运行命令
    try:
        run_command(aucell_cmd, logger)
    except subprocess.CalledProcessError as e:
        error_msg = f"pyscenic aucell command failed: {e}"
        if logger: logger.error(error_msg)
        else:
            print(error_msg)
        # 继续执行，因为即使没有aucell结果，我们也已经有了修剪后的网络
        aucell_output = None
    
    if aucell_output and os.path.exists(aucell_output):
        log_info(f"Regulon activity matrix saved to {aucell_output}")
        if logger: logger.debug(f"Regulon activity file size: {os.path.getsize(aucell_output)/1024/1024:.2f} MB")
    
    return ctx_output, aucell_output

# 分析结果
def analyze_results(pruned_grn_file, aucell_file, output_dir, output_prefix, logger=None):
    """分析GRN结果"""
    # 如果没有提供logger，使用默认的print函数
    if logger is None:
        log_info = print
    else:
        log_info = logger.info
    
    log_info("\n===== Step 3: Analyze GRN Results ======")
    
    # 加载修剪后的GRN
    log_info(f"Loading pruned GRN from {pruned_grn_file}...")
    if logger: 
        if os.path.exists(pruned_grn_file):
            logger.debug(f"Pruned GRN file size: {os.path.getsize(pruned_grn_file)/1024/1024:.2f} MB")
    pruned_grn = pd.read_csv(pruned_grn_file)
    if logger: logger.debug(f"Pruned GRN contains {len(pruned_grn)} interactions")
    
    # 统计调控子信息
    log_info("Analyzing regulon statistics...")
    regulon_counts = pruned_grn['regulon'].value_counts()
    log_info(f"Number of unique regulons: {len(regulon_counts)}")
    if logger: 
        logger.info(f"Top 10 regulons by gene count:")
        # 只记录前10个调控子到日志中
        for i, (regulon, count) in enumerate(regulon_counts.head(10).items()):
            logger.info(f"  {i+1}. {regulon}: {count} genes")
    else:
        print(f"Top 10 regulons by gene count:")
        print(regulon_counts.head(10))
    
    # 保存调控子统计信息
    regulon_stats_output = os.path.join(output_dir, f"{output_prefix}regulon_stats.txt")
    with open(regulon_stats_output, 'w') as f:
        f.write("Regulon Statistics:\n")
        f.write(f"Total regulons: {len(regulon_counts)}\n\n")
        f.write("Top 20 regulons by gene count:\n")
        for regulon, count in regulon_counts.head(20).items():
            f.write(f"{regulon}: {count} genes\n")
    
    log_info(f"Regulon statistics saved to {regulon_stats_output}")
    if logger: logger.debug(f"Regulon statistics file created")
    
    # 如果有aucell结果，分析调控子活性
    if aucell_file and os.path.exists(aucell_file):
        log_info(f"Loading regulon activity from {aucell_file}...")
        if logger: logger.debug(f"Regulon activity file size: {os.path.getsize(aucell_file)/1024/1024:.2f} MB")
        aucell_data = pd.read_csv(aucell_file, index_col=0)
        if logger: logger.debug(f"Regulon activity matrix shape: {aucell_data.shape}")
        
        # 计算每个调控子的平均活性
        regulon_activity = aucell_data.mean(axis=0)
        
        # 保存调控子活性统计
        activity_stats_output = os.path.join(output_dir, f"{output_prefix}regulon_activity_stats.txt")
        regulon_activity.sort_values(ascending=False).to_csv(activity_stats_output)
        
        log_info(f"Regulon activity statistics saved to {activity_stats_output}")
        if logger: 
            logger.debug(f"Regulon activity statistics file created")
            # 记录前5个最高活性的调控子
            top_5_active = regulon_activity.sort_values(ascending=False).head(5)
            logger.info(f"Top 5 most active regulons:")
            for i, (regulon, activity) in enumerate(top_5_active.items()):
                logger.info(f"  {i+1}. {regulon}: {activity:.4f}")
    
    # 生成总结报告
    report_output = os.path.join(output_dir, f"{output_prefix}analysis_report.txt")
    with open(report_output, 'w') as f:
        f.write("GRN Analysis Scheme 1 (scGLUE-based Pruning) Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"1. GRN Construction Details:\n")
        f.write(f"   - Pruned GRN file: {pruned_grn_file}\n")
        f.write(f"   - Number of regulons: {len(regulon_counts)}\n")
        f.write(f"   - Total gene-regulatory interactions: {len(pruned_grn)}\n\n")
        f.write(f"2. Key Findings:\n")
        f.write(f"   - Most frequent regulon: {regulon_counts.index[0]} ({regulon_counts.iloc[0]} genes)\n")
        f.write(f"   - Average genes per regulon: {regulon_counts.mean():.2f}\n\n")
        f.write(f"3. Files Generated:\n")
        f.write(f"   - Pruned GRN: {pruned_grn_file}\n")
        if aucell_file and os.path.exists(aucell_file):
            f.write(f"   - Regulon activity matrix: {aucell_file}\n")
        f.write(f"   - Regulon statistics: {regulon_stats_output}\n")
        if aucell_file and os.path.exists(aucell_file):
            f.write(f"   - Regulon activity statistics: {activity_stats_output}\n")
    
    log_info(f"Analysis report saved to {report_output}")
    if logger: logger.debug(f"Analysis report file created")

# 主函数
def main():
    args = parse_arguments()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 创建当前目录下的log目录
    log_dir = os.path.join(os.getcwd(), "log")
    os.makedirs(log_dir, exist_ok=True)
    
    # 创建日志文件
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"{args.output_prefix}log_{timestamp}.txt")
    
    # 设置日志，将日志级别设置为DEBUG以记录更详细信息
    logger = setup_logging(log_file, log_level=logging.DEBUG)
    logger.info("Starting GRN Analysis Scheme 1 process")
    logger.info(f"Command line arguments: {vars(args)}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Integration results directory: {args.integration_dir}")
    logger.info(f"Log file: {log_file}")
    logger.debug(f"Current working directory: {os.getcwd()}")
    logger.debug(f"Python version: {platform.python_version()}")
    
    try:
        # 步骤1：生成GRN草稿
        logger.info("===== Step 1: Generate GRN Draft ======")
        loom_file, grn_file = generate_grn_draft(
            args.rna_h5ad, args.tf_file, args.output_dir, args.output_prefix,
            args.num_workers, args.top_genes, logger
        )
        
        # 步骤2：使用scGLUE推断的feather文件修剪网络
        pruned_grn_file, aucell_file = prune_network_with_glue(
            grn_file, loom_file, args.integration_dir, args.output_dir, args.output_prefix,
            args.num_workers, args.rank_threshold, args.min_genes, args.feather_dir, logger
        )
        
        # 步骤3：分析结果
        analyze_results(pruned_grn_file, aucell_file, args.output_dir, args.output_prefix, logger)
        
        logger.info("\n===== GRN Analysis Scheme 1 Completed Successfully ======")
        logger.info(f"All results saved in: {args.output_dir}")
        logger.debug(f"Final file count in output directory: {len(os.listdir(args.output_dir))}")
    except Exception as e:
        logger.error(f"An error occurred during GRN Analysis Scheme 1: {str(e)}", exc_info=True)
        raise
    finally:
        logger.info("GRN Analysis Scheme 1 process completed")
        # 关闭所有日志处理器
        for handler in logging.getLogger().handlers[:]:
            handler.close()
            logging.getLogger().removeHandler(handler)

if __name__ == "__main__":
    main()