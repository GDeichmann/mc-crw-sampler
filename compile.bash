
mkdir obj
cd obj/
gcc ../src/*.c -I ../include/ -I $GMXLDLIB/../include/gromacs/ -c -g 
cd ../
gcc obj/*.o -L $GMXLDLIB -lm -lgmx -o mc_crw -g 

