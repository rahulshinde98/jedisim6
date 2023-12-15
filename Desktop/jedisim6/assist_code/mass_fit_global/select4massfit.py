# Read-in csv catalog
# Add column for z_source
# Add columns for g1 and g2
# Save new csv catalog



from astropy.table import Table
import numpy as np
import sys



#=======================
if len(sys.argv) != 3:
    print("python this.py csv_in_name csv_out_name")
    sys.exit(1)

#-----------------------


csv_in_name = sys.argv[1]
csv_out_name = sys.argv[2]
t = Table.read(csv_in_name, 
                format="ascii.csv")


#-----------------------
t["z"] = 0.6 * np.ones(len(t))
t["g1"] = t["e1"] / 2.0
t["g2"] = t["e2"] / 2.0


#-----------------------
print("x min, max: ", np.min(t["x"]), np.max(t["x"]))
print("y min, max: ", np.min(t["y"]), np.max(t["y"]))


#-----------------------
select = t["rmag"] > 14
select &= t["rmag"] < 26

select &= t["e1"] < 2
select &= t["e1"] > -2
select &= t["e2"] < 2
select &= t["e2"] > -2

t = t[select]

#-----------------------
print(t)
t.write(csv_out_name, 
        format="ascii.csv", overwrite=True)
