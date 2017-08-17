using QHull

p = rand(10,2)
ch = chull(p)
ch.points         # original points
ch.vertices       # indices to line segments forming the convex hull
ch.simplices      # the simplexes forming the convex hull
show(ch)
