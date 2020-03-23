dynmat.x < input_dynmat1.in > output_dynmat1
dynmat.x < input_dynmat2.in > output_dynmat2
dynmat.x < input_dynmat3.in > output_dynmat3
dynmat.x < input_dynmat4.in > output_dynmat4

cat dynmat1.eig dynmat2.eig dynmat3.eig dynmat4.eig > dynmat_all.eig
