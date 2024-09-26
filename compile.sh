# Bash files in ./dftb_scr/external/*.c

#
# Compile hessian.c file 
gcc -shared -fPIC -O3 hessian.c -o hessian.so

#
# Compile angular.c file
gcc -shared -fPIC -O3 angular.c -o angular.so


