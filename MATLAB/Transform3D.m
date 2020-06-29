function [rr] = Transform3D(r,X,Y,invert)
% Transform3D transforms a vector to a new set of orthonormal bases by
% calculating Tait-Bryan angles and performing the corresponding rotations
% Reference: https://en.wikipedia.org/wiki/Euler_angles
% Inputs:   r, vector in original basis frame
%           X, first basis vector
%           Y, second basis vector
%           invert, boolean value indicating whether
% Note:     Z, the third basis vector, is implied by X and Y
% Outputs:  rr, vector in new basis frame

if nargin < 4
    invert = false;
end

% normalize basis vectors
X = X/norm(X);
Y = Y/norm(Y);

psi = atan2(X(2),X(1));
theta = atan2(-X(3),sqrt(1-X(3)^2));
phi = atan2(Y(3),sqrt(1-Y(3)^2));

c1 = cos(psi);      s1 = sin(psi);
c2 = cos(theta);    s2 = sin(theta);
c3 = cos(phi);      s3 = sin(phi);

ZYX = [ c1*c2	c1*s2*s3-c3*s1	s1*s3+c1*c3*s2	;
        c2*s1   c1*c3+s1*s2*s3  c3*s1*s2-c1*s3  ;
        -s2     c2*s3           c2*c3           ];

if invert
    rr = ZYX*r;
else
    rr = ZYX'*r;
end

end
