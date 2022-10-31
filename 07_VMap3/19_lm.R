library(ggplot2)
a <- read.table("/Users/guoyafei/Desktop/test.txt",header=T,stringsAsFactors = F,sep="\t")

p <- ggplot(data = a, aes(x = Date, y = Proportion)) +
  geom_smooth(method = "lm",
   color="#558ebd",
                fill="lightgray",
               alpha=.7,
                size=0.8,se=T,
                formula = y~x) +
 geom_point(size=4,
            alpha=0.7,
           color="#fe654c")+ xlim(-10000,-6000)+ ylim(0,1)+
 theme_bw()
p

lm_eqn <- function(a){
  m <- lm(Proportion ~ Date, a);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 5),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

p1 <- p + 
  geom_text(x = -9000, y = 0.75, 
            label = lm_eqn(a), parse = TRUE)
p1
