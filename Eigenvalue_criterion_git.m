syms lambda11 lambda12 mu11 mu12 M1 M2 eps positive

a11 = -2/3*(3*lambda11*M1 + eps*lambda12*M2);
a12 = 2*(mu11+ (eps*lambda12*M1)/3);
a13 = -(2*eps*lambda12*M1)/(9*(1+mu12));
a14 = -(4*eps*lambda12)/(9*(1+mu12));
a21 = (lambda12*M2)/3 + (lambda11*M1)/eps; 
a22 = -mu11/eps-lambda12*M1/3;
a23 = -1/(2*eps)-(2*lambda12*M1)/(9*(1+mu12));
a24 = -4*lambda12/(9*(1+mu12));
a31 = 2* lambda11*M1/eps;
a32 = 0;
a33 = -(1/eps)*(2+mu11+((2*eps*lambda12*M1)*(3+mu12))/(3*(1+mu12)));
a34 = (2*lambda12*mu12)/(3*(1+mu12));
a41 = (2  * lambda12 * M1 * M2) / (9 * (1/2 + mu12)) + 2/eps * lambda11 * M1^2; 
a42 = 0;
a43 = -( lambda12 * M1^2) / (9 * (1 + mu12)) - 2/eps * M1;
a44 = -( ( lambda12 * M1 * (10 + 3 * mu12)) / (9 * (1 + mu12)) + mu11/eps );

A4 = 1;
A3 = - (a11 + a22 + a33 + a44); 
A2 = a11*a22 - a12*a21 + a11*a33 - a13*a31 + a11*a44 - a14*a41 + ...
     a22*a33 - a23*a32 + a22*a44 - a24*a42 + a33*a44 - a34*a43;
A1 = a11*a23*a32 - a11*a22*a33 + a12*a21*a33 - a12*a23*a31 - ...
     a13*a21*a32 + a13*a22*a31 - a11*a22*a44 + a11*a24*a42 + ...
     a12*a21*a44 - a12*a24*a41 - a14*a21*a42 + a14*a22*a41 - ...
     a11*a33*a44 + a11*a34*a43 + a13*a31*a44 - a13*a34*a41 - ...
     a14*a31*a43 + a14*a33*a41 - a22*a33*a44 + a22*a34*a43 + ...
     a23*a32*a44 - a23*a34*a42 - a24*a32*a43 + a24*a33*a42;
A0 = a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + ...
     a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - ...
     a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - ...
     a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + ...
     a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + ...
     a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - ...
     a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - ...
     a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41;

A2 = simplify(A2);
A1 = simplify(A1);
A0 = simplify(A0);

coeffs = [A4, A3, A2, A1, A0];

routh_array = sym([]);

routh_array(1, :) = [coeffs(1), coeffs(3), coeffs(5)];
routh_array(2, :) = [coeffs(2), coeffs(4), 0];

for i = 3:5
    for j = 1:2
        routh_array(i, j) = -( routh_array(i-2, 1)*routh_array(i-1, j+1) - ...
                               routh_array(i-2, j+1)*routh_array(i-1, 1) ) ...
                             / routh_array(i-1, 1);
    end
end

disp('Routh Array:');
first_col = simplify(routh_array(:, 1));

sign1 = sign(simplify(expand(A2*A3 - A1)));
D1 = simplify(expand(A1*A2*A3 - A1^2 - A0*A3^2));

assume(eps < 1)
assumeAlso(eps > 0)
assumeAlso(M1 > 0)
assumeAlso(M2 > 0)
assumeAlso(lambda11 > 0)
assumeAlso(lambda12 > 0)
assumeAlso(mu11 > 0)
assumeAlso(mu12 > 0)

sign2 = sign(D1);
sign3 = sign(first_col(5));
