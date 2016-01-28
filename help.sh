

KW=$1

grep -n -A 8 "$KW" read_in.f

 for BL in `grep -n -A 10 "$KW" read_in.f | grep -i "(S)" | awk '{print $2}' | sed "s/(S).*//g" `
 do
   echo "$BL" 
   echo ""
   grep "$BL" default.f
   echo ""
done

