Srcdfull.mmlcr1 <- mmlcr(outer = ~1 | id, components = list(
	list(formula = anti ~ poly(age, 2), min = 0, max = 12, 
	class = "cnormlong"), 
	list(formula = read ~ poly(age, 2), class = "normlong")), 
	data = Srcdfull, n.groups = 1)


Srcdfull.mmlcr2 <- update(Srcdfull.mmlcr1, n.groups = 2)
	
Srcdfull.mmlcr3 <- update(Srcdfull.mmlcr1, n.groups = 3)
	
Srcdfull.mmlcr4 <- update(Srcdfull.mmlcr1, n.groups = 4)
	
Srcdfull.mmlcr5 <- update(Srcdfull.mmlcr1, n.groups = 5)
	
	
anova(Srcdfull.mmlcr1, Srcdfull.mmlcr2, Srcdfull.mmlcr3, 
	Srcdfull.mmlcr4, Srcdfull.mmlcr5, test = "BIC")

plot(Srcdfull.mmlcr4, which = 1, smooth = 1)

plot(Srcdfull.mmlcr4, which = 2, smooth = 1)


Srcdfull.mmlcr4.homecog <- update(Srcdfull.mmlcr4, outer = ~homecog | id, 
	post.prob = post.prob(Srcdfull.mmlcr4))


Srcdfull.mmlcr4.homecogemo <- update(Srcdfull.mmlcr4, 
	outer = ~homecog + homeemo | id, 
	 post.prob = post.prob(Srcdfull.mmlcr4))
	
anova(Srcdfull.mmlcr4, Srcdfull.mmlcr4.homecog, Srcdfull.mmlcr4.homecogemo)

plot(fitted(Srcdfull.mmlcr4.homecogemo)[[1]], residuals(Srcdfull.mmlcr4.homecogemo)[[1]])
plot(fitted(Srcdfull.mmlcr4.homecogemo)[[2]], residuals(Srcdfull.mmlcr4.homecogemo)[[2]])