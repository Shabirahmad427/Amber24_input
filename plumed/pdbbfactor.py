with open("output.pdb", "r") as f_in, open("output.new.pdb", "w") as f_out:
    for line in f_in:
        if line.startswith("ATOM"):
            # Parse original residue number
            residue_num_str = line[22:26]
            try:
                residue_num = int(residue_num_str)
            except ValueError:
                residue_num = None

            # Change residue number 410 to 430
            if residue_num == 410:
                residue_num = 430
                new_res_num_str = f"{residue_num:>4}"
                line = line[:22] + new_res_num_str + line[26:]

            atom_name = line[12:16].strip()
            b_factor = line[60:66].strip()

            # Check if residue is in 175-190 or 410-430 and CA atom and B-factor is 0.00
            if (atom_name == "CA" and
                (175 <= residue_num <= 190 or 410 <= residue_num <= 430) and
                b_factor == "0.00"):
                line = line[:60] + "  1.00" + line[66:]

        f_out.write(line)
