#!/bin/bash
declare filename="largon.x"
declare -a files=("routines.f90"           \
                  "largon.f90" \
                  )

echo "##### Compilation log #####"
for file in *.x *.out; do
   if [ -f $file ]; then
      rm -f $file
      echo "||  --> $file removed"
   #else
   #   echo "  $file... nothing to do"
   fi
done

for i in "${files[@]}"
do
      gfortran -c "$i" 
      if [ $? -ne 0 ]
      then
          echo -e "|X| Compilation \e[31;5;82mfailed\e[0m for $i failed"
          exit 1
      else
          echo -e "|| Compilation \e[38;5;82msucceeded\e[0m for $i"
      fi
done

gfortran -fopenmp  "${files[@]}" -o "$filename"

for file in *.o ; do
   if [ -f $file ]; then
      rm -f $file
      echo "||  --> $file removed"
   #else
   #   echo "  $file... nothing to do"
   fi
done
echo "##### End of compilation #####"
