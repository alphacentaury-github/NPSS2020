      subroutine gauss(a,b,npoint,xri,wri)
      implicit real*8(a-h,o-z)
      real*8 xg(200),wg(200),xri(200),wri(200)
      call setmgl(npoint,xg,wg)
      do 20 j=1,npoint
      xri(j) = (a+b)/2.d0 + (b-a)/2.d0*xg(j)
      wri(j) = (b-a)/2.d0*wg(j)
   20 continue
      return
      end
      subroutine setmgl( n, points, weight )
      implicit real*8  ( a-h, o-z )
      real*8 points(200), weight(200)
      real*8 poin16(300)
      real*8, parameter :: pi = 3.14159265358979323846
      if (n.gt.200)  write (1,50)
50    format(' setmlg call with too many points')
      m = ( n + 1 ) / 2
      e1 = n * ( n + 1 )
      do 1 i = 1, m
      t = ( 4*i - 1 ) * pi / ( 4*n + 2 )
      x0 = ( 1.d0 - ( 1.d0 - 1.d0/n ) / ( 8.d0*n*n ) ) * cos(t)
      pkm1 = 1.d0
      pk = x0
      do 3 k = 2, n
      t1 = x0 * pk
      pkp1 = t1 - pkm1 - ( t1-pkm1 )/k + t1
      pkm1 = pk
      pk = pkp1
3     continue
      den = 1.d0 - x0*x0
      d1 = n * ( pkm1 - x0*pk )
      dpn = d1 / den
      d2pn = ( 2.d0*x0*dpn - e1*pk ) / den
      d3pn = ( 4.d0*x0*d2pn + (2.d0-e1)*dpn ) / den
      d4pn = ( 6.d0*x0*d3pn + (6.d0-e1)*d2pn ) / den
      u = pk / dpn
      v = d2pn / dpn
      h = -u * ( 1.d0 + 0.5d0*u*(v+u*(v*v-u*d3pn/(3.d0*dpn))))
      p = pk + h*(dpn+0.5d0*h*(d2pn+h/3.d0*(d3pn+0.25d0*h*d4pn)))
      dp = dpn + h*(d2pn+0.5d0*h*(d3pn+h*d4pn/3.d0))
      h = h - p / dp
      poin16(i) = x0 + h
      fx = d1 - h*e1*(pk+0.5d0*h*(dpn+h/3.d0*   
     &     (d2pn+0.25d0*h*(d3pn+0.2d0*h*d4pn))))  
      weight(i) = 2.d0 * ( 1.d0 - poin16(i)*poin16(i)) / (fx*fx)
1     continue
      if ( m + m .gt. n ) poin16(m) = 0.d0
      do 10 i = n/2 + 1, n
      poin16(i) = poin16( n + 1 - i )
      weight(i) = weight( n + 1 - i )
      poin16( n + 1 - i ) = -poin16( n + 1 - i )
10    continue
      do 30 i=1,n
 30   points(i)=poin16(i)
      return
      end
