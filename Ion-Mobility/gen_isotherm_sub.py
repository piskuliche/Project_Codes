#!/usr/bin/env python

def sbatch_lines(fname,nions):
    f = open(fname,'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --job-name=%d_ions\n' % nions)
    f.write('#SBATCH --partition=laird,thompson\n')
    f.write('#SBATCH --constraint=intel\n')
    f.write('#SBATCH --output=%d_ions_%%A_%%a.log\n'  % nions)
    f.write('#SBATCH --mail-user=piskuliche@ku.edu\n')
    f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH --mem=100G\n')
    f.write('#SBATCH --ntasks=20\n')
    f.write('#SBATCH --time=30:00:00\n')
    f.write('#SBATCH --array 0-10\n')
    f.close()
    return

def write_header(fname,nions):
    f = open(fname,'a')
    f.write('# Number of Ions\n')
    f.write('NCAT=%d\n' % nions)
    f.write('NANI=%d\n' % nions)
    f.write('\n')
    f.write('module load codecol/build\n')
    f.write('\n')
    f.write('mkdir -p $SLURM_ARRAY_TASK_ID/$NCAT\n')
    f.write('\n')
    f.close()
    return

def set_vars(fname, nacat, naani):
    f = open(fname, 'a')
    f.write('# Set the variables\n')
    f.write('LIN=$((SLURM_ARRAY_TASK_ID + 2))\n')
    f.write('X=$( awk \'NR == \'\"$LIN\"\' {print $1}\' yao.dat)\n')
    f.write('NCO2=$( awk \'NR == \'\"$LIN\"\' {print $2}\' yao.dat)\n')
    f.write('NACN=$( awk \'NR == \'\"$LIN\"\' {print $3}\' yao.dat)\n')
    f.write('LNG=$( awk \'NR == \'\"$LIN"\' {print $4}\' yao.dat)\n')
    f.write('\n')
    f.write('ACNSRTSKIP=9\n')
    f.write('ACNENDSKIP=$((NCO2*3 + NCAT*%d + NANI*%d))\n' % (nacat,naani))
    f.write('CO2SRTSKIP=$((9+NACN*3))\n')
    f.write('CO2ENDSKIP=$((NCAT*%d + NANI*%d))\n' % (nacat, naani))
    f.write('CATSTRTSKIP=$((NACN*3 + NCO2*3 + 9))\n')
    f.write('CATEND_SKIP=$((NANI*%d))\n' % naani)
    f.write('ANISRTSKIP=$((9 + NACN*3 + NCO2*3+NCAT*%d))\n' % nacat)
    f.write('ANIENDSKIP=0\n')
    f.close()
    return

def molinp_file(name, NAME):
    line = []
    line.append('echo "&nml" > %s.inp' % name)
    line.append('echo "nt=400," >> %s.inp' % name)
    line.append('echo "or_int=20," >> %s.inp' % name)
    line.append('echo "dt=0.1," >> %s.inp' % name)
    line.append('echo "startconfig=0," >> %s.inp' % name)
    line.append('echo "endconfig=500000," >> %s.inp' % name)
    line.append('echo "startskip=$%sSRTSKIP," >> %s.inp' % (NAME,name))
    line.append('echo "endskip=$%sENDSKIP," >> %s.inp' % (NAME,name))
    line.append('echo "nblocks=5," >> %s.inp' % name)
    line.append('echo "nmols=$N%s," >> %s.inp' % (NAME,name))
    line.append('echo "mol_name=\"acn\"," >> %s.inp' % name)
    line.append('echo "L=$LNG $LNG $LNG," >> %s.inp' % name)
    line.append('echo "needblock=.false." >> %s.inp' % name)
    line.append('echo "/" >> %s.inp' % name)
    return line

def build_directory(fname, mols):
    f = open(fname, 'a')
    f.write('\n')
    f.write('cd $SLURM_ARRAY_TASK_ID/$NCAT\n')
    f.write('\n')
    f.write('cp ../../build.inp ./\n')
    f.write('sed -i -e "s@AAA@$NACN@g" build.inp\n')
    f.write('sed -i -e "s@BBB@$NCO2@g" build.inp\n')
    f.write('sed -i -e "s@CCC@$NCAT@g" build.inp\n')
    f.write('sed -i -e "s@DDD@$NANI@g" build.inp\n')
    f.write('sed -i -e "s@EEE@$LNG@g" build.inp\n')
    f.write('cp ../../in.acncx ./\n')
    for i in range(4):
        f.write('cp ../../%s.txt ./\n' % mols[i])
    f.write('\n')
    f.write('build.py < build.inp\n')
    f.write('\n')
    f.write('cat lmps.paircoeffs lmps.bondcoeffs lmps.anglecoeffs > lmps.include.coeffs\n')
    f.write('\n')
    
    f.close()
    return

def write_lammps(fname):
    f = open(fname, 'a')
    f.write('module load lammps/.test12Dec2018\n')
    f.write('mpirun lmp_mpi -sf opt -in in.acncx\n')
    f.write('\n')
    f.close()

def write_msd(fname,mols,moltypes):
    f = open(fname, 'a')
    f.write('\n')
    f.write('# MSD Section\n')
    f.write('module load compiler/pgi/18\n')
    f.write('module load codecol/msd\n')
    f.write('export OMP_NUM_THREADS=20\n')
    f.write('export OMP_STACKSIZE=3G\n')
    f.write('rm msd_*.dat\n')
    f.write('rm cm_msd_*.dat\n')
    f.write('LNG=$(cat "NL.avg")\n')
    f.write('\n')
    for i in range(4):
        lines = molinp_file(mols[i], moltypes[i])
        for line in lines:
            f.write('%s\n' % line)
        f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('for i in {0..4}; do\n')
    f.write('\techo $i\n')
    f.write('\tcalc_msd.exe $i < acn.inp\n')
    f.write('\tcalc_msd.exe $i < co2.inp\n')
    f.write('\tcalc_msd.exe $i < pclo4.inp\n')
    f.write('\tcalc_msd.exe $i < li.inp\n')
    f.write('done\n')
    f.write('\n')
    f.write('line=$((SLURM_ARRAY_TASK_ID+1))\n')
    f.write('molfrac=$(sed -n "$line"p ../../xco2)\n')
    f.write('\n')
    for i in range(4):  
        f.write('parse_msds.py -f msd_%s -n 5 -startskip 100 -endskip 0 -molname %s -label $molfrac\n' % (mols[i],mols[i]))
    f.close()
    return

def write_atomfiles(mols):
    for i in range(4):
        f = open(mols[i]+'.txt', 'w')
        if mols[i]=='acn':
            natoms=3
            M=[15.035, 12.011, 14.007]
        elif mols[i]=='co2':
            natoms=3
            M=[15.999,12.011, 15.999]
        elif mols[i]=='li':
            natoms=1
            M=[6.94]
        elif mols[i]=='br':
            natoms=1
            M=[79.904]
        elif mols[i]=='cl':
            natoms=1
            M=[35.45]
        elif mols[i]=='pclo4':
            natoms=5
            M=[35.45,15.999,15.999,15.999,15.999]
        f.write('%d\n' % natoms)
        f.write('1\n')
        f.write('1 1\n')
        for m in M:
            f.write('%.4f\n' % m)
        f.close()

def write_ff(ff):
    f = open('build.inp','w')
    f.write('{\n')
    f.write('"num_components":4,\n')
    f.write('"nspec":[AAA,BBB,CCC,DDD],\n')
    f.write('"tspec":["%s", "%s", "%s", "%s"],\n' % (ff[0],ff[1],ff[2],ff[3]))
    f.write('"blength":EEE,\n')
    f.write('"e_unit":["kj","kcal","kcal","kcal"],\n')
    f.write('"f_unit":["kcal","kcal","kcal","kcal"],\n')
    f.write('"eo_unit":["kcal","kcal","kcal","kcal"],\n')
    f.write('"fo_unit":["kcal","kcal","kcal","kcal"],\n')
    f.write('"shift_f":1.0,\n')
    f.write('"ff_type":["lj","lj","lj","lj"]\n')
    f.write('}\n')
    f.close()

    


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', default='sub_NIONS.sh', help='Isotherm submit script name')
    parser.add_argument('-n', default=1, help='Number of ions')
    parser.add_argument('-m', default=[], help='Name of ions (should be 2)', action='append')
    parser.add_argument('-ff', default=[], help='Name of ff (should be 2)', action='append')
    args = parser.parse_args()


    nions = int(args.n)
    fname = str(args.f)
    mols = ['acn', 'co2']+args.m
    ff = ['acn_maroncelli', 'co2_epm2']+args.ff
    moltypes = ['ACN', 'CO2', 'CAT', 'ANI']
    write_ff(ff)
    sbatch_lines(fname, nions)
    write_header(fname,nions)
    set_vars(fname, 1,5)
    build_directory(fname, mols)
    write_msd(fname,mols,moltypes)
    write_atomfiles(mols)


