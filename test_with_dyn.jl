dawsonmia(x)=dawson(Float64(x))
expmia(x)=exp(big(x))
erfimia(x)=2/sqrt(pi)*expmia(big(x^2))*dawson(Float64(x))




x, w = gausslegendre(25)
t1, w1 = gausshermite(25)
t3, w3 = gausshermite(25)
global m1=0.1;
global q1=0.2;
global chi1=0.3;
global m3=0.1;
global q3=0.2;
global chi3=0.3;

function updateJD(j1,j3,B,bb,M)
	#global m1=rand();
	#global q1=rand();
	#global chi1=rand();
	#global m3=rand();
	#global q3=rand();
	#global chi3=rand();
	global h=1.0;
	global mu=M
	global beta=B;
	global betaprima=bb;
	global n=betaprima/beta;
	global J1=j1;
	global s1=1.0;
	global J1aster=beta*J1;
	global s1aster=beta*s1;
	global J3=j3;
	global s3=1/sqrt(betaprima*mu);
	#global s3=bb;
	global J3aster=beta*J3;
	global s3aster=beta*s3;
end

function update_param()
	c=s1aster/sqrt(2)*sqrt(abs(chi1-q1));
	d=s3aster/sqrt(2)*sqrt(abs(chi3-q3));;
	a(t1,t3)=-beta*h+J1aster*m1+s1aster*sqrt(q1)*t1-J3aster*m3-s3aster*sqrt(q3)*t3;
	b(t3)=beta*h-J3aster*m3-s3aster*sqrt(q3)*t3;
end

function converge(x,y,precision)::Bool
       if norm(x.-y)<precision
             return true
             else return false
       end
end


function converge_adv(x,y,z,precision)::Bool
       if converge(x,y,precision)==true
             return true
			 (m1,q1,chi1,m3,q3,chi3)=  (cumul1.+cumul2)./2
	   elseif converge(x,z,precision)==true
			 return true

	 (m1,q1,chi1,m3,q3,chi3)=  (cumul1.+cumul2)./2
	   else return false
	   end
end

global cumul1=zeros(6)
global cumul2=zeros(6)
global cumul3=zeros(6)


function calculations()
		  global m1_=m1;
		  global q1_=q1;
		  global chi1_=chi1;
		  global m3_=m3;
		  global q3_=q3;
		  global chi3_=chi3;
	global c=s1aster/sqrt(2)*sqrt(abs(chi1-q1));
	global d=s3aster/sqrt(2)*sqrt(abs(chi3-q3));
	 a(t1,t3)=-beta*h+J1aster*m1+s1aster*sqrt(q1)*t1-J3aster*m3-s3aster*sqrt(q3)*t3;
	 b(t3)=beta*h-J3aster*m3-s3aster*sqrt(q3)*t3;

	integrandZ(x,t1,t3)=(expmia((1/2*x-1/2)*(a(t1,t3)+(c^2+d^2)*(1/2*x-1/2)))*(-dawsonmia(b(t3)/(2*d)+d*(1/2*x-1/2))+expmia(b(t3)+d^2+2*d^2*(1/2*x-1/2))*dawsonmia(b(t3)/(2*d)+d+d*(1/2*x-1/2))))/d
	Z(t1,t3)=1/2*dot(w,integrandZ.(x,t1,t3))
	Z_n(t1,t3)=(Z(t1,t3))^n

	integrand_down(t1,t3)=Z_n(t1,t3)
	intermediate_down(t3)=1/sqrt(pi)*dot(w1,integrand_down.(sqrt(2)*t1,t3))
	global denom=1/sqrt(pi)*dot(w3,intermediate_down.(sqrt(2)*t3))

	#Fisrt_moment_1

	integrandZ_v1(x,t1,t3)=(expmia((1/2*x-1/2) *(a(t1,t3)+(c^2+d^2)*(1/2*x-1/2)))*(-dawsonmia(b(t3)/(2*d)+d*(1/2*x-1/2))+expmia(b(t3)+d^2+2*d^2*(1/2*x-1/2))*dawsonmia(b(t3)/(2*d)+d+d*(1/2*x-1/2))))*(1/2*x-1/2)/d
	Z_v1(t1,t3)=1/2*dot(w,integrandZ_v1.(x,t1,t3))

	integrand_up_v1(t1,t3)=Z_n(t1,t3)*Z_v1(t1,t3)/Z(t1,t3)
	intermediate_up_v1(t3)=1/sqrt(pi)*dot(w1,integrand_up_v1.(sqrt(2)*t1,t3))
	nominator_m_v1=1/sqrt(pi)*dot(w3,intermediate_up_v1.(sqrt(2)*t3))

	#Fisrt_moment_1_2

	integrand_up_v1_2(t1,t3)=Z_n(t1,t3)*(Z_v1(t1,t3)/Z(t1,t3))^2
	intermediate_up_v1_2(t3)=1/sqrt(pi)*dot(w1,integrand_up_v1_2.(sqrt(2)*t1,t3))
	nominator_q_v1=1/sqrt(pi)*dot(w3,intermediate_up_v1_2.(sqrt(2)*t3))

	#Second_moment

	integrandZ_v1square(x,t1,t3)=(expmia((1/2*x-1/2)*(a(t1,t3)+(c^2+d^2)*(1/2*x-1/2)))*(-dawsonmia(b(t3)/(2*d)+d*(1/2*x-1/2))+expmia(b(t3)+d^2+2*d^2*(1/2*x-1/2))*dawsonmia(b(t3)/(2*d)+d+d*(1/2*x-1/2))))*(1/2*x-1/2)^2/d
	Z_v1_square(t1,t3)=1/2*dot(w,integrandZ_v1square.(x,t1,t3))

	integrand_up_v1square(t1,t3)=Z_n(t1,t3)*Z_v1_square(t1,t3)/Z(t1,t3)
	intermediate_up_v1square(t3)=1/sqrt(pi)*dot(w1,integrand_up_v1square.(sqrt(2)*t1,t3))
	nominator_chi_v1=1/sqrt(pi)*dot(w3,intermediate_up_v1square.(sqrt(2)*t3))



	#Fisrt_moment_1_third

	integrandZ_v3(x,t1,t3)=-(expmia((1/2*x-1/2)*(a(t1,t3)+(c^2+d^2)*(1/2*x-1/2)))*(-d+b(t3)*dawsonmia(b(t3)/(2*d)+d*(1/2*x-1/2))+expmia(b(t3)+d^2+2*d^2*(1/2*x-1/2))*(d-b(t3)*dawsonmia(b(t3)/(2*d)+d+d*(1/2*x-1/2)))))/(2*d^3)
	Z_v3(t1,t3)=1/2*dot(w,integrandZ_v3.(x,t1,t3))
	integrand_up_v3(t1,t3)=Z_n(t1,t3)*Z_v3(t1,t3)/Z(t1,t3)
	intermediate_up_v3(t3)=1/sqrt(pi)*dot(w1,integrand_up_v3.(sqrt(2)*t1,t3))
	nominator_m_v3=1/sqrt(pi)*dot(w3,intermediate_up_v3.(sqrt(2)*t3))

	#Fisrt_moment_1_2

	integrand_up_v3_2(t1,t3)=Z_n(t1,t3)*(Z_v3(t1,t3)/Z(t1,t3))^2
	intermediate_up_v3_2(t3)=1/sqrt(pi)*dot(w1,integrand_up_v3_2.(sqrt(2)*t1,t3))
	nominator_q_v3=1/sqrt(pi)*dot(w3,intermediate_up_v3_2.(sqrt(2)*t3))

	#Second_moment

	integrandZ_v3square(x,t1,t3)=(1/(8*d^5))*expmia(-(b(t3)^2/(4*d^2))-b(t3)*(1/2*x-1/2)+(1/2*x-1/2)*(a(t1,t3)+c^2*(1/2*x-1/2)))*(2*d*expmia((b(t3)+2*d^2*(1/2*x-1/2))^2/(4*d^2))*(b(t3)-2*d^2*(1/2*x-1/2)+expmia(b(t3)+d^2+2*d^2*(1/2*x-1/2))*(-b(t3)+2*d^2*(1+(1/2*x-1/2))))-(b(t3)^2-2*d^2)*sqrt(pi)*(erfimia(b(t3)/(2*d)+d*(1/2*x-1/2))-erfimia(b(t3)/(2*d)+d+d*(1/2*x-1/2))))
	Z_v3_square(t1,t3)=1/2*dot(w,integrandZ_v3square.(x,t1,t3))

	integrand_up_v3square(t1,t3)=Z_n(t1,t3)*Z_v3_square(t1,t3)/Z(t1,t3)
	intermediate_up_v3square(t3)=1/sqrt(pi)*dot(w1,integrand_up_v3square.(sqrt(2)*t1,t3))
	nominator_chi_v3=1/sqrt(pi)*dot(w3,intermediate_up_v3square.(sqrt(2)*t3))

	global m10=nominator_m_v1/denom
	global q10=nominator_q_v1/denom
	global chi10=nominator_chi_v1/denom

	global m30=nominator_m_v3/denom
	global q30=nominator_q_v3/denom
	global chi30=nominator_chi_v3/denom

	if mod(counter,3)==0
		global cumul1=(m10,q10,chi10,m30,q30,chi30)
	elseif mod(counter,3)==1
		global cumul2=(m10,q10,chi10,m30,q30,chi30)
	else
		global cumul3=(m10,q10,chi10,m30,q30,chi30)
	end

end

trouble_par=[]
global datafull=[]
	global m1=-0.1;
	global q1=0.2;
	global chi1=0.3;
	global m3=0.1;
	global q3=0.2;
	global chi3=0.3;
function fixsolve(k,l,B,bb,mu)
	updateJD(k,l,B,bb,mu)
      global counter=0;
      test=false;
				#println((m10,q10,chi10))
			#	println((m30,q30,chi30))
            while test==false&&counter<100

				  calculations()
                  		  counter =counter +1

				  println("")
				  println(counter)
				  test=converge((m10,q10,chi10,m30,q30,chi30),(m1,q1,chi1,m3,q3,chi3),10.0^(-7.0))


				   global m1=(m10+m1)/2;
   				   global q1=(q10+q1)/2;
   				   global chi1=(chi10+chi1)/2;
   				   global m3=(m30+m3)/2;
   				   global q3=(q30+q3)/2;
   				   global chi3=(chi30+chi3)/2;

				   	   #global m1=(m1_+m1)/2;
					   #global q1=(q1_+q1)/2;
					   #global chi1=(chi1_+chi1)/2;
					   #global m3=(m3_+m3)/2;
					   #global q3=(q3_+q3)/2;
					   #global chi3=(chi3_+chi3)/2;

				  println((m1,q1,chi1))
				  println((m3,q3,chi3))

            end
            push!(datafull,[betaprima,beta,n,mu,h,J1,J3,s1,s3,m1,q1,chi1,m3,q3,chi3])
            writedlm(string("J3-8J10.dat"),datafull)
			#writedlm("trouble.dat",trouble_par)

end

infostring="betaprima,beta,n,mu,h,J1,J3,s1,s3,m1,q1,chi1,m3,q3,chi3"
writedlm("info.dat",infostring)

lista=collect(range(0.1, step=1.0, length=11)).^2
for i in 1:11
	lista[i]=1/lista[i]
end

for k in 0.001:0.0005:0.006
	for l in 1.0:0.5:8.0	
		println("mu=",k);
	     	println("n=",l);
			fixsolve(0.0,-10.0,0.05,0.05*l,k)
	end
end
#for i in 0.1:1.0:10
	#fixsolve(0,0,1.0,8.0)
#end
#for k in lista[1:11] 	 los k son el J1 original, mientras que J debajo es J1*
#	println("J1=",k);
#	println("s2_1=",1/sqrt(l));
#	for l in lista[1:11] #los l son el J3 original, mientras que J3 es J3*
#    	println("J3=",k);
#	println("s2_3=",1/sqrt(l));
#	fixsolve(k,l,1.0,2.0)
#	end
#end
