#This routine allows to modify params.py from the terminal

from aart_func import *

def nullable_string(val):
    if not val:
        return "None"
    return val

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='AART')

  parser.add_argument('--a', help='BH Spin e.g., 0.94')

  parser.add_argument('--i', help="Observers' inclination [degrees] e.g., 17")

  parser.add_argument('--maxbl', help="Max Baseline")

  parser.add_argument('--gamma', help="Astro param gamma")
  parser.add_argument('--mu', help="Astro param mu")
  parser.add_argument('--sigma', help="Astro param sigma")

  parser.add_argument('--subkep', help="Sub-Keplerianity")

  parser.add_argument('--betar', help="Radial velocity factor")

  parser.add_argument('--betaphi', help="Angular velocity factor")

  parser.add_argument('--snapshots', help="Number of Snapshots")

  parser.add_argument('--p_image', help="Make grids coincide")

  parser.add_argument('--limits', help="Limits of the image")

  parser.add_argument('--dx0', help="Resolution of n0")

  parser.add_argument('--dx1', help="Resolution of n1")

  parser.add_argument('--dx2', help="Resolution of n2")

  args = parser.parse_args()

  for line in fileinput.input("params.py", inplace=True):

    if line.strip().startswith('spin_case=') and args.a!=None:
        line = 'spin_case=%s\n'%args.a

    if line.strip().startswith('i_case=') and args.i!=None:
        line = 'i_case=%s\n'%args.i

    if line.strip().startswith('gammap=') and args.gamma!=None:
        line = 'gammap=%s\n'%args.gamma
    if line.strip().startswith('mup=') and args.mu!=None:
        line = 'mup=%s\n'%args.mu
    if line.strip().startswith('sigmap=')and args.sigma!=None:
      line = 'sigmap=%s\n'%args.sigma

    if line.strip().startswith('sub_kep=') and args.subkep!=None:
        line = 'sub_kep=%s\n'%args.subkep

    if line.strip().startswith('betar=') and args.betar!=None:
        line = 'betar=%s\n'%args.betar

    if line.strip().startswith('betaphi=') and args.betaphi!=None:
        line = 'betaphi=%s\n'%args.betaphi

    if line.strip().startswith('snapshots=') and args.snapshots!=None:
        line = 'snapshots=%s\n'%args.snapshots

    if line.strip().startswith('p_image=') and args.p_image!=None:
        line = 'p_image=%s\n'%args.p_image

    if line.strip().startswith('limits=') and args.limits!=None:
        line = 'limits=%s\n'%args.limits

    if line.strip().startswith('dx0=') and args.dx0!=None:
        line = 'dx0=%s\n'%args.dx0

    if line.strip().startswith('dx1=') and args.dx1!=None:
        line = 'dx1=%s\n'%args.dx1

    if line.strip().startswith('dx2=') and args.dx2!=None:
        line = 'dx2=%s\n'%args.dx2

    sys.stdout.write(line)

#print("The parameters were changed in params.py\n")

  #os.system('radialintensity intensity.py')