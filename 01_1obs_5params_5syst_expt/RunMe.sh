# TRENTo-based Bayesian analysis demo
# https://github.com/dereksoeder/TrentoBayesianDemo
#
# Last updated 05/02/2022


# Command-line help

if [ "$#" -ne 0 ]
then
  cat 1>&2 << ___
TRENTo-based Bayesian analysis demo
https://github.com/dereksoeder/TrentoBayesianDemo

This script will generate TRENTo events, process them into files usable by the
analysis code, and then run the analysis.  The most important output will be
'plots/posterior.pdf'.

Other than the TRENTo executable, no input is needed: just run '$0'.
___
  exit 1
fi


# Initialization

nevents=1000  # number of events per design point per collision system

trentopath='../trento/build/src/trento'
designfile='processed/Design.dat'

if ! [ -e "$trentopath" ]
then
  echo "'$trentopath' does not exist.  Please build TRENTo." 1>&2
  exit 1
fi

if ! [ -e "$designfile" ]
then
  echo "'$designfile' does not exist.  You might need to re-download or generate the design." 1>&2
  exit 1
fi

rm -f {,raw}badpoints.txt

numdespts=0

eventspath='data'
exptpath='expt'

if ! [ -e "$exptpath" ]
then
  echo "'$exptpath' does not exist.  The experimental data cannot be processed." 1>&2
  exit 1
fi

mkdir -p "$eventspath"


# For each system...

touch rawbadpoints.txt

for sys in 'Au Au 62d4' 'Au Au 200' 'Cu Cu 62d4' 'Cu Cu 200' 'Pb Pb 2760'
do
  read projectile target sqrts <<< "$sys"
  systag="${projectile}${target}${sqrts}"

  despt=0

  # For each design point...

  while read norm62d4 norm200 norm2760 nucleon_width fluctuation EXCESS
  do
    if ! ( [[ "${norm62d4}" =~ ^[0-9.Ee+-]+$ ]] && [[ "${norm200}" =~ ^[0-9.Ee+-]+$ ]] && [[ "${norm2760}" =~ ^[0-9.Ee+-]+$ ]] && [[ "${nucleon_width}" =~ ^[0-9.Ee+-]+$ ]] && [[ "${fluctuation}" =~ ^[0-9.Ee+-]+$ ]] && [ "$EXCESS" == '' ] )
    then
      continue  # skip unparseable lines such as header and blank lines
    fi

    # see http://qcd.phy.duke.edu/trento/usage.html
    if [ "$sqrts" == '62d4' ]
    then
      cross_section=3.60  # 1509.06727
      norm="${norm62d4}"
    elif [ "$sqrts" == '200' ]
    then
      cross_section=4.23  # 1509.06727
      norm="${norm200}"
    elif [ "$sqrts" == '2760' ]
    then
      cross_section=6.4  # 1108.6027
      norm="${norm2760}"
    else
      echo "Unsupported sqrt(s_NN) '$sqrts'." 1>&2
      exit 1
    fi

    # Generate events

    eventsfile="$eventspath/trento-${systag}-$(printf %03d $despt).txt"

    echo "Generating '$eventsfile'..." 1>&2

    "$trentopath" \
        "$projectile" "$target" "$nevents" \
        -x "${cross_section}" \
        -n "$norm" \
        -w "${nucleon_width}" \
        -k "$fluctuation" \
        -p 0. -d 0. \
        --grid-max 15. --grid-step 0.1 \
      > "$eventsfile"

    if [ "$?" -ne 0 ] || ( ! [ -e "$eventspath" ] ) || [ "$(cat "$eventsfile" | wc -l)" -ne "$nevents" ]
    then
      echo "WARNING: Event generation failed for design point ${despt}; marking it bad." 1>&2
      echo "$despt" >> rawbadpoints.txt
    fi

    despt="$(($despt + 1))"
  done < "$designfile"

  if [ "$numdespts" -lt "$despt" ]
  then
    numdespts="$despt"
  fi
done

if [ "$numdespts" -eq 0 ]
then
  echo "Event generation failed; '$designfile' might be invalid." 1>&2
  exit 1
fi

sort -g rawbadpoints.txt | uniq > badpoints.txt
rm rawbadpoints.txt


# Process events and experimental data

rm -f processed/{Data,Prediction}-*.dat

echo 1>&2
echo 'Processing data...' 1>&2

python3 Process.py

if [ "$?" -ne 0 ] || [ "$(readlink -e processed/Data-*.dat | wc -l)" -eq 0 ] || [ "$(readlink -e processed/Prediction-main-*.dat | wc -l)" -eq 0 ]
then
  echo "Process.py failed; analysis cannot continue.  Check that '$designfile' and 'badpoints.txt' are consistent and '$exptpath/*.dat' files are valid." 1>&2
  exit 1
fi


# Run Bayesian analysis

echo 1>&2
echo 'Running analysis...' 1>&2

# analysis code generates .hdf files in a 'cache' directory, and HDF5 by default tries to lock them
# either cache must reside on a file system that supports locking (on NERSC, scratch does but home does not), or locking must be disabled as below

HDF5_USE_FILE_LOCKING=FALSE \
python3 Analysis.py AuAu62d4,AuAu200,CuCu62d4,CuCu200,PbPb2760 ETmid 5 1500 500 1500
