using Plots 
function foo()
  pyplot()
  x = LinRange(0, 1, 100)
  y = exp.(-x.^2)
  p = plot(x, y)
  savefig(p, "media/test.pdf")
  display(p)
end