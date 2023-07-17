function weierstrass_form_cubic(f)
    x,y = gens(f.parent)
    t0 = coeff(f,y^3)
    s1 = coeff(f,y^2)
    s0 = coeff(f,y^2*x)
    r2 = coeff(f,y)
    r1 = coeff(f,x*y)
    r0 = coeff(f,y*x^2)
    q3 = coeff(f, f.parent(1))
    q2 = coeff(f,x)
    q1 = coeff(f,x^2)
    q0 = coeff(f,x^3)
    if (t0==0) && (s1==1) && (s0==0) && (r0==0) && (q0==-1)
        return f 
    end
    a1=r1;
    a2=-(s0*q2+s1*q1+r0*r2);
    a3=(9*t0*q0-s0*r0)*q3+((-t0*q1-s1*r0)*q2+(-s0*r2*q1-s1*r2*q0));
    a4=((-3*t0*r0+s0^2)*q1+(-3*s1*s0*q0+s1*r0^2))*q3+(t0*r0*q2^2+(s1*s0*q1+((-3*t0*r2+s1^2)*q0+s0*r0*r2))*q2+(t0*r2*q1^2+s1*r0*r2*q1+s0*r2^2*q0)); 
    a6=(-27*t0^2*q0^2+(9*t0*s0*r0-s0^3)*q0-t0*r0^3)*q3^2+(((9*t0^2*q0-t0*s0*r0)*q1+((-3*t0*s0*r1+(3*t0*s1*r0+2*s1*s0^2))*q0+(t0*r0^2*r1-s1*s0*r0^2)))*q2+(-t0^2*q1^3
    +(t0*s0*r1+(2*t0*s1*r0-s1*s0^2))*q1^2+((3*t0*s0*r2+
	(-3*t0*s1*r1+2*s1^2*s0))*q0+((2*t0*r0^2-s0^2*r0)*r2+
	(-t0*r0*r1^2+s1*s0*r0*r1-s1^2*r0^2)))*q1+((9*t0*s1*r2-
	s1^3)*q0^2+(((-3*t0*r0+s0^2)*r1-s1*s0*r0)*r2+(t0*r1^3
	-s1*s0*r1^2+s1^2*r0*r1))*q0)))*q3+(-t0^2*q0*q2^3+
	(-t0*s1*r0*q1+((2*t0*s0*r2+(t0*s1*r1-s1^2*s0))*q0-
	t0*r0^2*r2))*q2^2+(-t0*s0*r2*q1^2+(-t0*s1*r2*q0+
	(t0*r0*r1-s1*s0*r0)*r2)*q1+((2*t0*r0-s0^2)*r2^2+
	(-t0*r1^2+s1*s0*r1-s1^2*r0)*r2)*q0)*q2+
	(-t0*r0*r2^2*q1^2+(t0*r1-s1*s0)*r2^2*q0*q1-
     t0*r2^3*q0^2));
    b2=a1^2+4*a2;
    b4=2*a4+a1*a3;
    b6=a3^2+4*a6;
    b8=a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;
    c4=b2^2-24*b4;
    c6=-b2^3+36*b2*b4-216*b6;
    return(y^2+a1*x*y+a3*y-x^3-a2*x^2-a4*x-a6)
end

function j_invariant_cubic(f)
  f=weierstrass_form_cubic(f)
  x,y = gens(f.parent)
  a1 = coeff(f,x*y)
  a2 = -coeff(f,x^2)
  a3 = coeff(f,y)
  a4 = -coeff(f,x)
  a6 = -coeff(f,f.parent(1))
  b2=a1^2+4*a2
  b4=2*a4+a1*a3
  b6=a3^2+4*a6
  b8=a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2
  c4=b2^2-24*b4
  delta=-b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6
  if (delta==0)
    print("The input is a rational curve and has no j-invariant!")
    return
  end
  invariant=(div(c4^3,delta))
  return(invariant)
end



