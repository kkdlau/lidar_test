#include "Circle.hpp"
#include "PointArray.hpp"
#include "Statistics.hpp"
#include "helper_funcs.hpp"

namespace CircleFittingModel {
	void eigen2x2(double a, double b, double c, double& d1, double& d2, double& Vx, double& Vy)
/*            <------- Input ----------> <--------------- Output ----------------->

       Eigendecomposition of a symmetric 2x2 matrix
          faster and more accurate than the library function
       
       Input:  a,b,c - components of the matrix [a c
                                                 c b]
       Output:  d1,d2 - eigenvalues
                Vx,Vy - unit eigenvector for d1
                
                The eigenvector for d2 need not be computed, it is (-Vy,Vx)
                
       Note:  d1 is the leading eigenvalue, i.e., |d1|>=|d2|
       
       Nikolai Chernov,  June 2012

*/
{
    double disc,f;

    disc = pythag(a-b,2.0*c);    // discriminant

    d1 = (a+b > 0.) ? (a + b + disc)/2.0 : (a + b - disc)/2.0;
    d2 = (a*b - c*c)/d1;

    if (abs(a-d1) > abs(b-d1))
    {
        if ((f=pythag(c,d1-a))==0.) 
        {
            Vx = One; Vy = 0.;  return;
        }
        else       
        {
            Vx = c/f;  Vy = (d1-a)/f;
        }
    }
    else
    {
        if ((f=pythag(c,d1-b))==0.) 
        {
            Vx = One; Vy = 0.;  return;
        }
        else       
        {
            Vx = (d1-b)/f;  Vy = c/f;
        }
    }
    
    return;
}

double SigmaWithLargeCircleOption (PointArray<double>& data, Circle& circle)
/*                                <-------- Input --------->

		Compute the objective function	for the geometric circle fitting problem
		
		Input:  data     - the class of data (contains the given points):
		
		        data.n   - the number of data points
		        data.X[] - the array of X-coordinates
		        data.Y[] - the array of Y-coordinates
		        data.meanX - the mean X-coordinate
		        data.meanY - the mean Y-coordinate
		          (the last 2.0 must be precomputed)
		          
		        circle    - the class of circle parameters:
		        
		        circle.a - the X-coordinate of the center
		        circle.b - the Y-coordinate of the center
		        
		Output:
		        the value of the objective function
		        (more precisely, the square root of the average square of the distance) 
          
		Nikolai Chernov,  January 2013
*/
{
    int i,n=data.n;
    double sum=0.,dx,dy,r,D[n];
    double LargeCircle=Ten,a0,b0,del,s,c,x,y,z,p,t,g,W,Z;
   
    if (abs(circle.a)<LargeCircle && abs(circle.b)<LargeCircle)   // main case (not a large circle)
    {
    	for (i=0; i<n; i++)
    	{
    		dx = data.X[i] - circle.a;
    		dy = data.Y[i] - circle.b;
    		D[i] = sqrt(dx*dx+dy*dy);
    		sum += D[i];
    	}
    	r = sum/n;
    	
    	for (sum=0., i=0; i<n; i++)  sum += SQR(D[i] - r);
    	
    	return sum/n;
    }
    else   //  case of a large circle
    {
    	a0 = circle.a-data.meanX;  b0 = circle.b-data.meanY;
    	del = One/sqrt(a0*a0 + b0*b0);
    	s = b0*del;  c = a0*del;
    	
    	for (W=Z=0.,i=0; i<n; i++)
    	{
    		x = data.X[i] - data.meanX;
    		y = data.Y[i] - data.meanY;
    		z = x*x + y*y;
    		p = x*c + y*s;
    		t = del*z - 2.0*p;
    		g = t/(One+sqrt(One+del*t));
    		W += (z+p*g)/(2.0+del*g);
    		Z += z;
    	}
    	W /= n;
    	Z /= n;
    	
    	return Z-W*(2.0+del*del*W);
    }
}

void GradientHessian (Data& data, Circle& circle, double& F1, double& F2, double& A11, double& A22, double& A12)
/*                    <-------- Input --------->  <----------------------- Output ----------------------->

		Compute the gradient vector and the Hessian matrix of the objective function
		        for the geometric circle fitting problem
		
		Input:  data     - the class of data (contains the given points):
		
		        data.n   - the number of data points
		        data.X[] - the array of X-coordinates
		        data.Y[] - the array of Y-coordinates
		        data.meanX - the mean X-coordinate
		        data.meanY - the mean Y-coordinate
		          (the last 2.0 must be precomputed)
		          
		        circle    - the class of circle parameters:
		        
		        circle.a - the X-coordinate of the center
		        circle.b - the Y-coordinate of the center
		        
		Output:
		        [F1 F2]  - the coordinates of the gradient vector
		        
		        A11 A12  - the components of the Hessian matrix
		        A12 A22    (it is symmetric, so only three are computed)
          
		Nikolai Chernov,  January 2013
*/
{
	int i,n=data.n;
    double LargeCircle=Ten,dx,dy,r,u,v,Mr,Mu,Mv,Muu,Mvv,Muv,Muur,Mvvr,Muvr;
	double a0,b0,del,dd,s,c,x,y,a,b,z,p,t,w,g,g1,gg1,gg2;
	double X,Y,R,U,V,T,W,AA,BB,AB,AG,BG,GG,UUR,VVR,UVR;
    
    if (abs(circle.a)<LargeCircle && abs(circle.b)<LargeCircle)   // main case (not a large circle)
    {
    	for (Mr=Mu=Mv=Muu=Mvv=Muv=Muur=Mvvr=Muvr=0.,i=0; i<n; i++)
    	{
    		dx = data.X[i] - circle.a;
    		dy = data.Y[i] - circle.b;
    		r = sqrt(dx*dx + dy*dy);
    		u = dx/r;
    		v = dy/r;
    		Mr += r;
    		Mu += u;
    		Mv += v;
    		Muu += u*u;
    		Mvv += v*v;
    		Muv += u*v;
    		Muur += u*u/r;
    		Mvvr += v*v/r;
    		Muvr += u*v/r;
    	}
    	Mr /= n;
    	Mu /= n;
    	Mv /= n;
    	Muu /= n;
    	Mvv /= n;
    	Muv /= n;
    	Muur /= n;
    	Mvvr /= n;
    	Muvr /= n;
    	
    	F1 = circle.a + Mu*Mr - data.meanX;
    	F2 = circle.b + Mv*Mr - data.meanY;
    	
    	A11 = One - Mu*Mu - Mr*Mvvr;
    	A22 = One - Mv*Mv - Mr*Muur;
    	A12 = -Mu*Mv + Mr*Muvr;
    }
    else   //  case of a large circle
    {
    	
    	a0 = circle.a-data.meanX;  b0 = circle.b-data.meanY;
    	del = One/sqrt(a0*a0 + b0*b0);  dd = del*del;
    	s = b0*del;  c = a0*del;
    	
    	for (X=Y=R=T=W=AA=BB=AB=AG=BG=GG=0.,i=0; i<n; i++)
    	{
    		x = data.X[i] - data.meanX;
    		y = data.Y[i] - data.meanY;
    		z = x*x + y*y;
    		p = x*c + y*s;
    		t = 2.0*p-del*z;
    		w = sqrt(One-del*t);
    		g = -t/(One+w);
    		g1 = One/(One+del*g);
    		gg1 = g*g1;
    		gg2 = g/(2.0+del*g);
    		a = (x+g*c)/w;
    		b = (y+g*s)/w;
    		X += x*gg1;
    		Y += y*gg1;
    		R += z + t*gg2;
    		T += t*gg1;
    		W += t*gg1*gg2;
    		AA += a*a*g1;
    		BB += b*b*g1;
    		AB += a*b*g1;
    		AG += a*gg1;
    		BG += b*gg1;
    		GG += g*gg1;
    	}
    	X /= n;
    	Y /= n;
    	R /= n;
    	T /= n;
    	W /= n;
    	AA /= n;
    	BB /= n;
    	AB /= n;
    	AG /= n;
    	BG /= n;
    	GG /= n;
    	
    	U = (T-del*W)*c/2.0 - X + R*c/2.0; 
    	V = (T-del*W)*s/2.0 - Y + R*s/2.0;
    	
//         compute the components of the gradient vector
    	
    	F1 = del*((dd*R*U - del*W*c + T*c)/2.0 - X);
    	F2 = del*((dd*R*V - del*W*s + T*s)/2.0 - Y);
    	
    	UUR = ((GG-R/2.0)*c + 2.0*(AG-U))*c + AA;
    	VVR = ((GG-R/2.0)*s + 2.0*(BG-V))*s + BB;
    	UVR = ((GG-R/2.0)*c + (AG-U))*s + (BG-V)*c + AB;
    	
//         compute the components of the Hessian matrix
    	
    	A11 = dd*(U*(2.0*c - dd*U) - R*s*s/2.0 - VVR*(One + dd*R/2.0));
    	A22 = dd*(V*(2.0*s - dd*V) - R*c*c/2.0 - UUR*(One + dd*R/2.0));
    	A12 = dd*(U*s + V*c + R*s*c/2.0 - dd*U*V + UVR*(One + dd*R/2.0));
	}
}

int CircleFitByChernovHoussam (Data& data, Circle& circleIni, double LambdaIni, Circle& circle)
/*                             <------------------ Input ------------------->  <-- Output -->

       Geometric circle fit to a given set of data points (in 2D)
		
       Input:  data     - the class of data (contains the given points):
		
	       data.n   - the number of data points
	       data.X[] - the array of X-coordinates
	       data.Y[] - the array of Y-coordinates
		          
               circleIni - parameters of the initial circle ("initial guess")
		        
	       circleIni.a - the X-coordinate of the center of the initial circle
	       circleIni.b - the Y-coordinate of the center of the initial circle
	       circleIni.r - the radius of the initial circle
		        
	       LambdaIni - the initial value of the control parameter "lambda"
	                   for the Levenberg-Marquardt procedure
	                   (common choice is a small positive number, e.g. 0.001)
		        
       Output:
	       integer function value is a code:
	                  0:  normal termination, the best fitting circle is 
	                      successfully found
	                  1:  the number of outer iterations exceeds the limit (99)
	                      (indicator of a possible divergence)
	                  2:  the number of inner iterations exceeds the limit (99)
	                      (another indicator of a possible divergence)
		          3:  convergence to a point where the Hessian matrix 
		              is NOT positive definite. This indicates that 
		              the fitting circle may correspond to a maximum
		              or to a saddle point of the objective function
		                   
	       circle - parameters of the fitting circle ("best fit")
		        
	       circle.a - the X-coordinate of the center of the fitting circle
	       circle.b - the Y-coordinate of the center of the fitting circle
 	       circle.r - the radius of the fitting circle
 	       circle.s - the root mean square error (the estimate of sigma)
 	       circle.i - the total number of outer iterations (updating the parameters)
 	       circle.j - the total number of inner iterations (adjusting lambda)
 	       
       Algorithm is based on the paper by H. Abdul-Rahman and N. Chernov
       "Fast and numerically stable circle fit" 
       in Journal of Mathematical Imaging and Vision (2013)
          
		Nikolai Chernov,  January 2013
*/
{
    int i,n=data.n,iter,inner,IterMAX=200,check_line=1,code;
    
    double lambda;
    double F1,F2,A11,A22,A12,dX,dY,Mxx,Myy,Mxy,Mxxy,dx,dy;
    double d1,d2,dmin=One,Vx,Vy,dL1,dL2,VLx,VLy,aL,bL,R,G1,G2,sBest,gBest,AB,progress;

//          control parameters (have been optimized empirically):

    double ParLimit2=100.;
    double epsilon=1.e+7*REAL_EPSILON;;
    double factor1=32.,factor2=32.;
    double ccc=0.4;
    double factorUp=10.,factorDown=0.1;
    
    Circle Old,New;
    
    data.means();   // Compute x- and y-means (via a function in class "data") 
    
    //    starting with the given initial guess
    
    New = circleIni;            //  initialize the circle
    New.s = SigmaWithLargeCircleOption(data,New);    //  compute the root-mean-square error
    GradientHessian(data,New,F1,F2,A11,A22,A12);  // compute the gradient vector and Hessian matrix
    New.Gx = F1;  New.Gy = F2;   New.g = sqrt(F1*F1 + F2*F2);   //  the gradient vector and its norm 
    
    lambda = LambdaIni;         //    initialize lambda
    iter = inner = code = 0;    //    initialize iteration counters and the exit code
    sBest = gBest = progress = REAL_MAX;   //  set control variables to LARGE values
    //if (lpr==1) cout << iter <<"  ("<<New.a<<","<<New.b<<")  s="<<New.s<<"  g="<< New.g<<"  L="<<lambda << endl;
    
NextIteration:   //  starting point of the current iteration of Newton's method
	
    if (iter>0) progress = (abs(New.a-Old.a)+abs(New.b-Old.b))/(SQR(Old.a)+SQR(Old.b)+One);
               // evaluate the progress made during the previous iteration
    Old = New;
    if (++iter > IterMAX) goto enough;   //  termination due to going over the limit

    eigen2x2(A11,A22,A12,d1,d2,Vx,Vy);  //  eigendecomposition of the Hessian matrix
    dmin =  (d1<d2) ? d1 : d2;          //  recording the smaller e-value
    
	AB=sqrt(SQR(Old.a)+SQR(Old.b)) + One;   //  approximation to the circle size
	
//     main stopping rule: terminate iterations if 
//          the gradient vector is small enough and the progress is not substantial 
    if ((Old.g < factor1*REAL_EPSILON)&&(progress<epsilon))
    {
    		//if (lpr==1) cout << "++++ gradient is small enough ++++" << endl;
    		goto enough;
    }
    
//     secondary stopping rule (prevents some stupid cycling)
    if ((Old.s >= sBest)&&(Old.g >= gBest))
    {
    		//if (lpr==1) cout << "++++ iterations stabilize (best results repeated) ++++" << endl;
    		goto enough;
    }
   
    if (sBest > Old.s) sBest = Old.s;  //  updating the smallest value of the objective function found so far
    if (gBest > Old.g) gBest = Old.g;  //  updating the smallest length of the gradient vector found so far
	
	G1 = Vx*F1 + Vy*F2;  //  rotating the gradient vector
	G2 = Vx*F2 - Vy*F1;  //  (expressing it in the eigensystem of the Hessian matrix)

try_again:   //  starting point of an "inner" iteration (adjusting lambda)

//           enforcing a lower bound on lambda that guarantees that
//           (i)  the augmented Hessian matrix is positive definite
//           (ii) the step is not too big (does not exceed a certain fraction of the circle size)
//                                         the fraction is defined by the factor "ccc")
	if (lambda < abs(G1)/AB/ccc - d1)  lambda = abs(G1)/AB/ccc - d1;
	if (lambda < abs(G2)/AB/ccc - d2)  lambda = abs(G2)/AB/ccc - d2;

//           computing the step (dX,dY) by using the current va;ue of lambda

    dX = Old.Gx*(Vx*Vx/(d1+lambda)+Vy*Vy/(d2+lambda)) + Old.Gy*Vx*Vy*(One/(d1+lambda)-One/(d2+lambda));
    dY = Old.Gx*Vx*Vy*(One/(d1+lambda)-One/(d2+lambda)) + Old.Gy*(Vx*Vx/(d2+lambda)+Vy*Vy/(d1+lambda));
   
//           updating the circle parameters

    New.a = Old.a - dX;
    New.b = Old.b - dY;
    
    if ((New.a==Old.a)&&(New.b==Old.b))   // if no change, terminate iterations    
    {
    		//if (lpr==1) cout << "++++ iterations stabilize (no change in center) ++++" << endl;
    		goto enough;
    }

//       check if the circle is very large

    if (abs(New.a)>ParLimit2 || abs(New.b)>ParLimit2)
    {
//          when the circle is very large for the first time, check if 
//          the best fitting line gives the best fit

    		if (check_line)  //  initially, check_line=1, then it is set to zero
    		{
    			//if (lpr==1) cout << "  Linear zone 1st  iter=" << iter << "  (" << New.a << "," << New.b << ")" << endl;
    		
   //                compute scatter matrix
   
  	  		for (Mxx=Myy=Mxy=0.,i=0; i<n; i++)   
    			{
    				dx = data.X[i] - data.meanX;
    				dy = data.Y[i] - data.meanY;
    				Mxx += dx*dx;
    				Myy += dy*dy;
    				Mxy += dx*dy;
    			}
    	
    			eigen2x2(Mxx,Myy,Mxy,dL1,dL2,VLx,VLy);  //  eigendecomposition of scatter matrix

//                   compute the third mixed moment (after rotation of coordinates)

    			for (Mxxy=0.,i=0; i<n; i++)
    			{
    				dx = (data.X[i] - data.meanX)*VLx + (data.Y[i] - data.meanY)*VLy;
    				dy = (data.Y[i] - data.meanY)*VLx - (data.X[i] - data.meanX)*VLy;
    				Mxxy += dx*dx*dy;
    			}
//              check if the line is the best fit

    			//if (Mxxy == 0.) { cout << "  Line " << endl; cin.ignore(); }  //  need to finish this block...

//              rough estimate of the center to be used later to recoved from the wrong valley

			R = (Mxxy>0.) ? ParLimit2 : -ParLimit2;
			aL = -VLy*R;
			bL =  VLx*R;                 
    			check_line = 0;              //  set to zero to prevent further checks for line 
    		}
		
		if ((New.a*VLy - New.b*VLx)*R>0.)  // check if the circle is in the wrong valley
		{
	    		//if (lpr==1) cout << "  Linear zone foul  iter=" << iter << "  (" << New.a << "," << New.b << ")" << endl;
			New.a = aL;                        //  switch to the rough circle estimate
			New.b = bL;                        //    (precomupted earlier)
			New.s = SigmaWithLargeCircleOption(data,New);           //  compute the root-mean-square error
			GradientHessian(data,New,F1,F2,A11,A22,A12);  // compute the gradient vector and Hessian matrix
			New.Gx = F1;  New.Gy = F2;   New.g = sqrt(F1*F1 + F2*F2);   //  the gradient vector and its norm 
	    		lambda = LambdaIni;                //  reset lambda
	    		sBest = gBest = REAL_MAX;          //  reset best circle characteristics 
	    		//if (lpr==1) cout << "  Linear zone flip  iter=" << iter << "  (" << New.a << "," << New.b << ")" << endl;
			goto NextIteration;      //  restart the Newton's iteration
		}
    }
    	
    New.s = SigmaWithLargeCircleOption(data,New);      //  compute the root-mean-square error
    GradientHessian(data,New,F1,F2,A11,A22,A12);  // compute the gradient vector and Hessian matrix
    New.Gx = F1;  New.Gy = F2;   New.g = sqrt(F1*F1 + F2*F2);   //  the gradient vector and its norm 
    
    //if (lpr==1) cout << setprecision(15)<<iter <<"  ("<<New.a<<","<<New.b<<"  s="<<New.s<<"  g="<< New.g<<"  L="<<lambda<<endl;
    	
//                check if improvement is gained
    	
    if (New.s < sBest*(One+factor2*REAL_EPSILON))    //    yes, improvement
    {
    		lambda *= factorDown;     //  reduce lambda
    		goto NextIteration;       //  proceed to the next Newton's iteration
    }
    else                                             //  no improvement
    {
    		//if (lpr==1) cout << "   repeat with higher lambda" << endl;
    		if (++inner > IterMAX) goto enough;   //  termination due to going over the limit
    		lambda = factorUp*lambda;             //  increace lambda
    		goto try_again;                       //  make another inner iteration                 
    }
    
enough:                   //           exit
	
    Old.r = OptimalRadius(data,Old);
    Old.i = iter;
    Old.j = inner;
    
    circle = Old;    //  output circle class
    
    if (iter  > IterMAX) code = 1;    //  error code set to 1
    if (inner > IterMAX) code = 2;    //  error code set to 2
    if ((dmin <= 0.)&&(code==0)) 
    { 
    		//cout << " negative e-value=" << setprecision(20) << dmin << " iter=" << iter <<"  ("<<Old.a<<","<<Old.b<< ")" <<  endl; 
    		code = 3;     //  error code set to 3
    }
    return code;
}



}