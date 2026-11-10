import math
def transform(matrix,point):
    if len(matrix)!=6:raise ValueError('仿射矩阵需要六个数')
    a,b,tx,c,d,ty=matrix;x,y=point
    return [a*x+b*y+tx,c*x+d*y+ty]
def compose(after,before):
    a,b,x,c,d,y=after;e,f,u,g,h,v=before
    return [a*e+b*g,a*f+b*h,a*u+b*v+x,c*e+d*g,c*f+d*h,c*u+d*v+y]
def signed_area(ring):
    if len(ring)<3:raise ValueError('环至少三个点')
    return sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(ring,ring[1:]+ring[:1]))/2
def segment_length(a,b):return math.hypot(b[0]-a[0],b[1]-a[1])
def interpolate(a,b,t):
    if not 0<=t<=1:raise ValueError('参数越界')
    return [x+(y-x)*t for x,y in zip(a,b)]
