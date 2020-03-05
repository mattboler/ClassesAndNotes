clc; clear all; close all;

%{
Matt Boler
HW3
%}

%% Constants
G = [1 1 1; 0 1 1; 0 0 1];
G_down = G' * G;
G_up = G * G';

%% 1. Consider G

%% a. What is the eigen decomposition of G?
[V, D] = eig(G)

% Is it significant?

% Assuming significance means nonzero, distinct eigenvalues, then no, as all eigenvalues are the same.

%% b. Show:
[V, Gamma, U] = svd(G)

% G = V * Gamma * U
disp("G = V * Gamma * U")
G
V * Gamma * U'

% G_down = G' * G
disp("G_down = G' * G")
G_down
G' * G

% V is the eigen-matrix of G_up
disp("V is the eigen-matrix of G_up")
V
[temp_v, temp_d] = eig(G_up);
temp_v

% U is the eigen-matrix of G_down
disp("U is the eigen-matrix of G_down")
U
[temp_v, temp_d] = eig(G_down);
temp_v

% Gamma's nonzero values are equal to the square root of the eigenvalues of
% G_up and G_down
disp("Gamma's nonzero values are equal to the square root of the eigenvalues of G_up and G_down")
eig_G_up = eig(G_up);
eig_G_down = eig(G_down);
Gamma
sqrt(eig_G_down)
sqrt(eig_G_up)

% G_up = G*G'
disp("G_up = G*G'")
G_up
G*G'

%% 2. What is null of a matrix if the matrix has distinct, non-zero eigenvalues?

% The null of the matrix would be just the 0 vector.

%% 3. Show null(G' . G) = null(G), R(G' . G) = R(G'), null(G.G') = null(G'), and R(G.G') = R(G)

%{
Important definitions:
colspace : range(A)
kernel : null(A)
coimage : range(A')
cokernel : null(A')

For a matrix A = V * Gamma * U' with rank r :
range(A) = first r columns of V
null(A) = last n-r columns of U'
range(A') = first r columns of U'
null(A') = last m-r columns of V
%}

% Important matrices:
G_prime_G = G' * G;
G_G_prime = G*G';

% null(G' . G) = null(G)
disp("null(G' . G) = null(G)")
[m1, n1] = size(G_prime_G);
r1 = rank(G_prime_G);

[m2, n2] = size(G);
r2 = rank(G);

[v1, g1, u1] = svd(G_prime_G);
u1 = u1';
[v2, v2, u2] = svd(G);
u2 = u2';

if (n1 - r1) == 0 && (n2 - r2) == 0
    disp('Rank of both null spaces is 0 -> same null space');
else
    disp("These span the same space")
    b1 = u1(:, (n1-r1):end)
    b2 = u2(:, (n2-r2):end)
end



% R(G' . G) = R(G')

disp("R(G' . G) = R(G')")

[m1, n1] = size(G_prime_G);
r1 = rank(G_prime_G);

[m2, n2] = size(G');
r2 = rank(G');

[v1, g1, u1] = svd(G_prime_G);
u1 = u1';
[v2, v2, u2] = svd(G');
u2 = u2';

if r1 == 0 && r2 == 0
    disp('Rank of both ranges is 0 -> same range space');
else
    disp("These span the same space")
    b1 = v1(:, 1:r1)
    b2 = v2(:, 1:r2)
end

% null(G.G') = null(G')
disp("null(G.G') = null(G')")
[m1, n1] = size(G_G_prime);
r1 = rank(G_G_prime);

[m2, n2] = size(G');
r2 = rank(G');

[v1, g1, u1] = svd(G_G_prime);
u1 = u1';
[v2, v2, u2] = svd(G');
u2 = u2';

if (n1 - r1) == 0 && (n2 - r2) == 0
    disp('Rank of both null spaces is 0 -> same null space');
else
    disp("These span the same space")
    b1 = u1(:, (n1-r1):end)
    b2 = u2(:, (n2-r2):end)
end

% R(G.G') = R(G)

disp("R(G.G') = R(G)")

[m1, n1] = size(G_G_prime);
r1 = rank(G_G_prime);

[m2, n2] = size(G);
r2 = rank(G);

[v1, g1, u1] = svd(G_G_prime);
u1 = u1';
[v2, v2, u2] = svd(G);
u2 = u2';

if r1 == 0 && r2 == 0
    disp('Rank of both ranges is 0 -> same range space');
else
    disp("These span the same space")
    b1 = v1(:, 1:r1)
    b2 = v2(:, 1:r2)
end

%% 4. Show the minimum error inverse solution satisfies v* orthogonal to null(G)

% Since we solve for a vector which when multiplied by G gives us u.
% If v_bar is in the null space of G, then it will always give 0 instead of
% u. So v_bar cannot be in the null space of G.

