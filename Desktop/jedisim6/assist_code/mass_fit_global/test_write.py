import os
script_path = os.path.abspath(__file__)
script_dir = os.path.split(script_path)[0]
rel_path =  "output_files/number.txt"
num_file_path = os.path.join(script_dir, rel_path)
f1 = open(num_file_path,'r')
num = f1.readlines()[0]
outputs = "output_files/outputs" + num.split()[0] + ".txt"
abs_file_path = os.path.join(script_dir, outputs)
f2 = open(abs_file_path, 'w')
for i in range(10):
    f2.write('%e\n'%i)
f1.close()
f2.close()
