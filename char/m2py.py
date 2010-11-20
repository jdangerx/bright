import numpy as np
import re


if_idx_str = ("""if (exist("idx", "var"));\n"""
              """  idx = idx + 1;\n"""
              """else;\n"""
              """  idx = 1;\n"""
              """end;"""
              )

num_pattern = "([0-9]+[.]?[0-9]*[Ee]?[+-]?[0-9]*)"

numpy_array_pattern = r"\[[0-9\sEe+-.,%]*\]"
matlab_array_pattern = r"\[[0-9\sEe+-.%]*\]"

lhs_variable_pattern = r"(\w+)\s*(\(idx.*?\))"
rhs_variable_pattern = r"(\w+)\s*\(idx.*?\)\s*=\s*(.*)"

def convert_res(filename):
    """Convert a matlab *.m file to a python file."""
    with open(filename, 'r') as mfile:
        f = mfile.read()

    # Keep comments around
    f = f.replace('%', '#')

    # Grab the number of 'if' statements
    IDX = f.count(if_idx_str)

    # Replace if statements with something more meaningful
    fpart = f.partition(if_idx_str)
    f = fpart[0] + "idx = 0" + fpart[2]
    f = f.replace(if_idx_str, 'idx += 1')

    # Replace matlab arrays with python lists
    arrays = re.findall(matlab_array_pattern, f)
    for a in arrays:
        new_a = re.sub(num_pattern, lambda mo: mo.group(0) + ',', a)
        f = f.replace(a, new_a)

    # Encapsulate python lists in numpy arrays
    f = re.sub(numpy_array_pattern, lambda mo: 'np.array(' + mo.group(0) + ')', f)


    # Add imports to header
    header = "import numpy as np\n\n"

    # Find all variables and shape
    vars_shape = np.unique( re.findall(lhs_variable_pattern, f) )
    vars_dtype = dict( re.findall(rhs_variable_pattern, f) )
    # Initialize variables to zero
    header = header + "# Initialize variables\n"
    for vs in vars_shape:
        # Determine shape
        s = re.search(r'\[.*:(.*?)\]', vs[1])
        if s == None:
            vs_shape = ""
        else:
            vs_shape = s.group(1)
            vs_shape = vs_shape.split()
            vs_shape = ", ".join(vs_shape)

        # Determine Data type
        rhs = vars_dtype[vs[0]]
        if ("\'" in rhs) or ("\"" in rhs):
            dt = "'S{0}'".format(int( s.group(1) ))
            vs_shape = ""
        elif ('.' in rhs) or ('E' in rhs) or ('e' in rhs):
            dt = "float"
        else:
            dt = "int"

        zero_line = "{0} = np.zeros([{1}, {2}], dtype={3})\n".format(vs[0], IDX, vs_shape, dt)
        header = header + zero_line

    # Add IDx to file
    if 0 < IDX:
        header = header + "\n\n# Maximum Index\n\nIDX = {0}\n\n".format(IDX)

    # Add header to file
    f = header + f

    # Replace variable overrides
    vars = np.unique( re.findall("(" + lhs_variable_pattern + ")", f) )
    for v in vars:
        f = f.replace(v[0], "{0}[idx] ".format(v[1]))

    # Write the file out
    new_filename = filename.rpartition('.')[0] + '.py'
    with open(new_filename, 'w') as pyfile:
        pyfile.write(f)

    return f
