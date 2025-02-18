      subroutine band(a,m,n,p)
      integer g,h,i,j,k,m,n,p(n),r
      complex a(n,1),c
      real eps,max,d
      r = (m+1)/2
      eps = 1.0e-10
      do 10 i = 1,n
      do 11 j = m+1,m+r-1
      a(i,j) = cmplx(0.0,0.0)
   11 continue
   10 continue

      do 20 k = 1,n
c
c     Find pivots
c
      max = 0.0
      i = k
      j = r
   25 if(i.gt.n .or. j.lt.1) goto 30
      d = cabs(a(i,j))
      if(max.ge.d) goto 35
      max = d
      p(k) = i
   35 i = i+1
      j = j-1
      goto 25
   30 continue
      if(max.le.eps) goto 99
c
c     Switch pivot rows
c
      if (p(k).eq.k) goto 40
        i = r
        j = r+k-p(k)
   50   if (i.gt.m+r-1 .or. i.gt.n-k+r) goto 40
        c = a(k,i)
        a(k,i) = a(p(k),j)
        a(p(k),j) = c
        i = i+1
        j = j+1
        goto 50
   40 continue
c
c     Decompose A
c
      a(k,r) = 1.0/a(k,r)
      h = r-1
      i=k+1
   60 if(h.lt.1 .or. i.gt.n) goto 20
      a(i,h) = a(i,h) * a(k,r)
        j = h+1
        g = r+1
   70   if(g.gt.m+r-1 .or. j.gt.n+r-i) goto 80
        a(i,j) = a(i,j) - a(i,h)*a(k,g)
        j = j+1
        g = g+1
        goto 70
   80 continue
      i = i+1
      h = h-1
      goto 60
   20 continue
      return
99    m = 0
      return
      end