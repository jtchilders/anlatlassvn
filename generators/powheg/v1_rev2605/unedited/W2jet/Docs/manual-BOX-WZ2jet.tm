<TeXmacs|1.0.7.14>

<style|article>

<\body>
  <doc-data|<doc-title|The POWHEG-BOX-W2jet and Z2jet manual>||>

  <\with|par-mode|center>
    <with|font-series|light|J.M. Campbell, R.K. Ellis, P. Nason and G.
    Zanderighi>
  </with>

  <section|Introduction>

  The <with|font-family|tt|POWHEG-BOX-W2jet> and <with|font-family|tt|Z2jet>
  programs <cite|Campbell:2013vha> are NLO generators for the QCD production
  of <math|W<rsup|\<pm\>>jj> and <math|Z jj> events in hadronic collisions,
  with the vector boson decaying into leptons. They can be properly
  interfaced to Shower Monte Carlo programs. Off-shell effects are \ included
  in the calculation. The vector bosons are given a finite, fixed width. The
  CKM matrix is assumed diagonal.\ 

  \ \ This document describes the input parameters that are specific to this
  implementation, and also certain new features of the
  <with|font-family|tt|POWHEG-BOX> that they use. The parameters that are
  common to all <with|font-family|tt|POWHEG BOX> implementations are given in
  the <with|font-family|tt|manual-BOX.pdf> document, in the
  <with|font-family|tt|POWHEG-BOX/Docs> directory.

  <section|Running the program>

  Enter the process directory and do<next-line><with|font-family|tt|$ make
  pwhg_main><next-line>to build the executable for the generation of the
  events. In the following, the input parameters to be entered in the
  <with|font-family|tt|powheg.input> file are described.

  <section|Process specific input parameters>

  Parameters in <with|font-family|tt|powheg.input> that are specific to
  <math|W jj> production:<next-line><with|font-family|tt|vdecaymodeW -11
  \ \ \ ! decay mode of W+/W- (11=e-,13=mu-,15=tau-,
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ !
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -11=e+, -13=mu+,
  -15=tau+.)><next-line><with|font-family|tt|mllmin \ \ \ \ 0 \ \ \ \ \ \ !
  Minimum invariant mass for the lepton pair (default 0)<next-line>bwcutoff
  \ \ 15 \ \ \ \ \ ! How many W widths above and below its pole mass do
  we<next-line> \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! allow for the
  virtuality of the W (default: absent)><next-line><next-line>Parameters in
  <with|font-family|tt|powheg.input> that are specific to <math|Z jj>
  production:<next-line><with|font-family|tt|vdecaymodeZ 11 \ \ \ \ ! Z decay
  mode: 11=e,12=nue, etc.>,<next-line><with|font-family|tt|mllmin \ \ \ \ 1
  \ \ \ \ \ \ ! Min invariant mass for the lepton pair (default 1
  GeV)<next-line>mllmax \ \ \ \ 140 \ \ \ \ ! Max invariant mass for the
  lepton pair (default sqrt(S))<next-line>bwcutoff \ \ 15 \ \ \ \ \ ! How
  many Z widths above and below its pole mass do we<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! allow for the virtuality of the Z
  (default: absent)><next-line>We remark that the restrictions on the
  invariant mass of the lepton pairs are all applied
  simultaneously.<next-line><next-line>The following parameters are the same
  for <math|W jj> and <math|Z jj>:<next-line><with|font-family|tt|runningscales
  0 \ \ \ ! if /=1 or absent use MW, if 1 uses
  HT/2><next-line><with|font-family|tt|bornsuppfact \ 20 \ \ ! Set the pt
  scale parameter for Born-zero damping factor><next-line><with|font-family|tt|
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! equal to 20 GeV. If absent, no
  damping factor><line-break><next-line><with|font-family|tt|ptborncut \ 0.01
  \ \ \ ! pt of jets and their relative pt limited by ptborncut<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! default: 0.01 GeV><next-line>

  <section|Option that activate the new <with|font-family|tt|POWHEG BOX>
  features>

  The <with|font-family|tt|W2jet> and <with|font-family|tt|Z2jet> packages
  use a number of new features that will eventually become the default in an
  upcoming version 2 of the <with|font-family|tt|POWHEG BOX>. Several new
  files, stored in the directory <with|font-family|tt|POWHEG-BOX/Version-pre2-1>,
  replace the files in the <with|font-family|tt|POWHEG-BOX> directory. The
  <with|font-family|tt|Makefile> is structured in such a way that the new
  files have priority over the old ones. Here we describe these new features.

  <\itemize-minus>
    <item><with|font-family|tt|MiNLO>: \ the procedure of ref.
    <cite|Hamilton:2012np>. Using <with|font-family|tt|MiNLO> also inclusive
    quantities not requiring that the jets are resolved are computed with a
    certain accuracy.<cite*|><next-line><with|font-family|tt|minlo \ \ \ 1
    \ \ \ \ \ ! Activate the minlo feature (default absent)><next-line>

    <item><with|font-family|tt|doublefsr 1 \ \ \ \ ! Activate doubling of
    singular region in final state<next-line> \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ !
    quark gluon splitting. Default on for this process.><next-line>This
    feature is described in ref. <cite|Nason:2013uba>.

    <item><with|font-family|tt|olddij \ \ \ 1 \ \ \ \ ! Use old values for
    the d_ij functions (default 0)><next-line>The new separation or regions
    described in ref. <cite|Campbell:2013vha> is used by default. To fall
    back on the traditional approach, set this flag to 1.

    <item><with|font-family|tt|fastbtlbound \ \ \ 1 \ \ \ \ ! Use improved
    upper bounding envelope (default 0)><next-line>This uses the method
    described in ref. <cite|Melia:2011gk>.

    <item>Reweighting: it is possible to run a process saving reweighting
    information in the event file, and then use this information to compute
    new weights corresponding, for example, to a different choice of the
    factorization and renormalisation scales. In order for reweighting to
    work one should always turn on the option<next-line><with|font-family|tt|withnegweights
    1><next-line>(on by default in Z2jet). In order to store reweighting
    information, one should include the line:<next-line><with|font-family|tt|storeinfo_rwgt
    \ 1 ! Store reweighting information (default 1)><next-line>If this
    feature is present during the event generation, a line of reweighting
    information will appear at the end of each event in the event file, of
    the form<next-line><with|font-family|tt|#rwgt 1 59 759.386766596014 100
    1001 0><next-line>At the end of the run, one can change a few parameters
    in the <with|font-family|tt|powheg.input> file, or even alter the
    <with|font-family|tt|POWHEG> executable, in a way that affects the value
    of the cross section (typically one can change scale factors, PDF's,
    etc.), include the line<next-line><with|font-family|tt|compute_rwgt 1
    \ \ \ ! compute new weight><next-line><with|font-family|tt|>and run the
    program again in exactly the same directory that was used to produce the
    event files. Similar event files, with the <with|font-family|tt|-rwgt>
    string in the name, will be produced, where at the end of each event a
    new weight is appended, of the form<next-line><with|font-family|tt|#new
    weight,renfact,facfact,pdf1,pdf2 19335.9 1.0 1.0 10050 10050 lha
    ><next-line>where the new weight, the new values for the renormalization
    and factorization scale factors, the pdf number for each hadron, and the
    pdf package that was used are reported. The reweighting run is much
    faster than the generation run. Thus, iterating this method, it is easy
    to append a bunch of lines with new weights to the end of each event, in
    order to perform uncertainty studies. Notice that if the program is run
    with no changed parameters, the new weight should equal the original
    weight (the third entry of the first line in the Les Houches event
    record) up to truncation errors.
  </itemize-minus>

  Besides the features listed above, the method used to perform parallel runs
  has been completely overhauled and improved. In the following subsection we
  describe how it works.

  <subsection|Parallel runs>

  It is possible to split a <with|font-family|tt|POWHEG> run into several
  parallel ones. However, since <with|font-family|tt|POWHEG> runs in several
  stages, each stage has to be run until the end before the next stage is
  started, in order for the parallelization to be useful. We need the
  following flags:<next-line><with|font-family|tt|manyseeds 1 ! get the seeds
  for the random number generator from the file
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ !
  pwgseeds.dat.<next-line>parallelstage \ \<less\>m\<gtr\> \ !
  \<less\>m\<gtr\>=1...4, which level of parallel stage
  <next-line>xgriditeration \<less\>n\<gtr\> \ ! \<less\>n\<gtr\> is the
  iteration level for the calculation of the<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! importance sampling grid
  improvement (relevant only for<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! parallelstage=1)><next-line>The
  <with|font-family|tt|pwgseeds.dat> file contains a list of integer seeds,
  one per line. If the <with|font-family|tt|manyseeds> flag is set, when the
  <with|font-family|tt|pwhg_main> program starts, it asks for an integer
  number <math|m>. This number is used to select a line in the
  <with|font-family|tt|pwgseeds.dat> file, and the integer found there is
  used to initialize the random number generator. In this way, several runs
  can be performed simultaneously in the same directory, all using different
  random number seeds. Furthermore, the integer <math|m> appears as a four
  digit integer (with leading zeros) in the name of all files that are
  generated by the corresponding run.

  The <with|font-family|tt|pwhg_main> program runs in four stages:

  <\enumerate-numeric>
    <item>An importance sampling grid for the integration of the inclusive
    cross section is built and stored. This stage is run when
    <with|font-family|tt|parallelstage=1>. The number of calls to the
    inclusive cross section is controlled by the input variable
    <with|font-family|tt|ncall1>. If also a remnant cross section is present,
    the number of calls for the remnants is equal to
    <with|font-family|tt|ncall1>, unless a parameter
    <with|font-family|tt|ncall1rm> appears. Initially, this stage is called
    with <with|font-family|tt|xgriditeration=1>. The information that is
    needed to compute the importance sampling grids are stored in the files
    named (assuming <math|m=1>) <with|font-family|tt|pwggridinfo-btl-xg1-0001.dat>
    for the <math|<wide|B|~>> function and
    \ <with|font-family|tt|pwggridinfo-rmn-xg1-0001.dat> for the remnant. The
    quality of the grid for each single run is represented in the topdrawer
    files <with|font-family|tt|pwg-xg1-0001-btlgrid.top>. These plots are not
    very significant, since what counts is the grid formed by assembling all
    the information contained in all the <with|font-family|tt|pwggridinfo>
    files. These are generated in the subsequent step, are named
    <with|font-family|tt|pwg-0001-btlgrid.top>, and are all equal among each
    other for all values of <math|m> that have been used.

    Observe that it is important that all
    \ <with|font-family|tt|parallelstage=1> runs complete, in order to use
    all the produced grid information for building the grid.

    If after the <with|font-family|tt|parallelstage=1> and
    \ <with|font-family|tt|xgriditeration=1> run the grids do not look smooth
    enough, one or more iterations with <with|font-family|tt|xgriditeration=2>
    and higher, can be attempted. At this stage, the relevant files are named
    <with|font-family|tt|pwggridinfo-btl-xg2-0001.dat>,
    <with|font-family|tt|pwggridinfo-rmn-xg2-0001.dat> and
    <with|font-family|tt|pwg-xg2-0001-btlgrid.top>, and so on. Increasing
    \ <with|font-family|tt|ncall1> also helps in getting more satisfactory
    grids.

    <item>This stage is performed after the
    <with|font-family|tt|parallelstage=1> run, by setting
    <with|font-family|tt|parallelstage=2>. The integral of the inclusive
    cross section is computed, and an upper bounding of the integrand, having
    the form of a multidimensional step function on the grid, is computed and
    stored. The number of calls used for this task is controlled by <math|>
    <with|font-family|tt|ncall2> and <with|font-family|tt|itmx2> (the total
    number of calls per run is <with|font-family|tt|ncall2*itmx2>). Again,
    independent variables for the remnants are also available if needed,
    \ <with|font-family|tt|ncall2rm> and <with|font-family|tt|itmx2rm>. If
    the files with information on the importance sampling grid are missing,
    the program complains and stops. The integration and upper bounding
    envelope information is stored in files named
    <with|font-family|tt|pwggrid-0001.dat>. If the
    \ <with|font-family|tt|storemintupb> flag is set, auxiliary files named
    <with|font-family|tt|pwgbtildeupb-0001.dat> and
    <with|font-family|tt|pwgremnupb-0001.dat> are created and loaded with
    these information. Again, all processes run at this stage should be
    completed before going on. Only at the next stage these files are all
    loaded and assembled to get the cross section and upper bounding
    envelopes using the full statistics of the runs. Files named
    <with|font-family|tt|pwgfullgrid-0001.dat> are then created. They are all
    equal among each other. On subsequent runs, if any
    \ <with|font-family|tt|pwgfullgrid-0001.dat> is present, it is loaded in
    place of all the others, since it contains all the necessary information.

    <item>This is obtained by setting <with|font-family|tt|parallelstage=3>.
    The upper bounding factors for the generation of radiation are computed
    at this stage. The number of calls (per run) is controlled by the
    variable <with|font-family|tt|nubound>. <math|>The computed factors are
    stored in the files named \ \ <with|font-family|tt|pwgubound-0001.dat>.

    <item>With <with|font-family|tt|parallelstage=4>, the program starts
    generating events, that are stored in files named
    <with|font-family|tt|pwgevents-0001.lhe>. The number of events per file
    is controlled by the \ <with|font-family|tt|numevts> parameter in the
    <with|font-family|tt|powheg.input> file.
  </enumerate-numeric>

  When performing parallel runs, for very slow processes that require very
  many cpu's, the number of events per process may not be enough for building
  a valid upper bounding grid (see ref. <cite|Melia:2011gk>). In this case,
  it is convenient to set the flag \ <with|font-family|tt|storemintupb=1>. In
  this way, the results of the calls to the inclusive cross section are all
  stored in files, and the upper bounding envelope is built at a later
  parallel stage, by reading and assembling all these files.\ 

  <section|Files to inspect at the end of the run>

  The file <with|font-family|tt|pwgstat.dat> (for single run), or
  <with|font-family|tt|pwg-stY-0001-stat.dat>, with the stage number
  <with|font-family|tt|Y> equal to 2 or 3, contain information on the value
  of the cross section. Notice that, in case of parallel runs, the stage 3
  files are the final ones, since they refer to the result of the combination
  of all the stage 2 runs, while the stage 2 files only include the result of
  a single run. It is important to check that the cross section has the
  required precision.

  At the end of the runs, the values of certain counters are printed in files
  named <with|font-family|tt|pwgcounters.dat>, or
  <with|font-family|tt|pwgcounters-stY-XXXX.dat>, where
  <with|font-family|tt|Y> is the <with|font-family|tt|parallelstage> value
  and <with|font-family|tt|XXXX> stands for the four digits run number in
  case of parallel runs. It is important to check in the stage 4 files
  (<with|font-family|tt|pwgcounters-st4-XXXX.dat>) that the number of upper
  bound violation in the generation of radiation and in the inclusive cross
  section are much smaller than the total number of events. They appears as
  follows:<with|font-family|tt|<next-line>upper bound failures in generation
  of radiation = 636.000000000000<next-line>upper bound failure in inclusive
  cross section = 340.0000000000000><next-line>These numbers should be of the
  order of a percent of the total number of events. If the number of failure
  for the inclusive cross section is too large, the number of calls
  <with|font-family|tt|ncalls2> should be increased. If the failure in the
  generation of radiation is too large, then \ <with|font-family|tt|nubound>
  should be increased. In this last case, one can also attempt to increase
  \ <with|font-family|tt|xupbound>, with no need to rerun the stage of
  generation of the upper bounds for radiation (that is to say, stage 3).

  The counter files also report the approximate value of the total number of
  seconds spent doing the virtual and real calculations. These numbers are
  significant only if the program is run at 100% speed. If the virtual time
  is large compared to the real time, it may be convenient to increase the
  folding parameters <with|font-family|tt|foldcsi>,
  \ <with|font-family|tt|foldy> and \ <with|font-family|tt|foldphi>, that
  also reduces the fraction of negative weighted events (for the allowed
  values for these parameter see the <with|font-family|tt|POWHEG BOX>
  manual).

  <section|Shower and analysis>

  The user can shower the Les Houches event with whatever software he likes.
  We provide elementary interfaces to <with|font-family|tt|PYTHIA> and
  <with|font-family|tt|Pythia8>. An interface to <with|font-family|tt|HERWIG>
  can be easily developed on the basis of the examples in other
  <with|font-family|tt|POWHEG BOX> processes. An enhancement of the previous
  histogramming package in the <with|font-family|tt|POWHEG BOX> does now
  output the histograms associated with different weights, if present in the
  Les Houches event file. Assuming that all event files are merged into a
  single <with|font-family|tt|pwgevents.lhe> files, in order to run the
  builtin analysis one proceeds as follows:<next-line><with|font-family|tt|$
  cd POWHEG-BOX/W2jet><next-line><with|font-family|tt|$ make
  main-PYTHIA-lhef><next-line><with|font-family|tt|$ cd
  test-lhc><next-line><with|font-family|tt|$
  ../main-PYTHIA-lhef><next-line>where <with|font-family|tt|test-lhc> is the
  directory where the events are stored. The output files
  are<next-line><with|font-family|tt|pwgPOWHEG+PYTHIA-output.top.><next-line>If
  more weights are present in the event file, files of the
  form<next-line><with|font-family|tt|pwgPOWHEG+PYTHIA-output-WY.top><next-line>are
  produced, where <with|font-family|tt|Y> is the weight ordering number in
  the event file. These files are gnuplot data files, and can be plotted with
  the gnuplot program.

  The analysis can also be run in parallel. If the program does not find a
  \ <with|font-family|tt|pwgevents.lhe> file, it prompts for a filename. If
  the file name has the form <with|font-family|tt|pwgevents-XXXX.lhe>, then
  the output files have the name<next-line><with|font-family|tt|pwgPOWHEG+PYTHIA-output-XXXX-WY.top><next-line>This
  is useful if one wants to speed up the analysis of parallel runs. It is up
  to the user, at the end, to combine the files. A fortran program
  <with|font-family|tt|mergedata.f> is provided for this purpose in the
  <with|font-family|tt|Version-pre2-1> directory. Its use is self
  explanatory.

  <section|Examples>

  In the directories <with|font-family|tt|runs-paper> we have put the setup
  to perform parallel runs of the kind that we used for our paper
  <cite|Campbell:2013vha>. This setup assumes that the run is to be performed
  interactively on a multi-core machine with 48 cores. The script
  <with|font-family|tt|runpar.sh> executes all stages of the run. It is up to
  the user to adapt the run to his/her batch environment. The example of the
  <math|W jj> takes about 2 hours for the preparation stage, and 10 hours for
  the generation of events, on a 48 cores machine. For the <math|Z jj> case,
  the preparation stage takes roughly eight hours, while the events are
  generated in about 75 hours.

  <\bibliography|bib|JHEP|paper.bib>
    <\bib-list|1>
      <bibitem*|1><label|bib-Campbell:2013vha>J.<nbsp>M. Campbell, R.<nbsp>K.
      Ellis, P.<nbsp>Nason, and G.<nbsp>Zanderighi, <with|font-shape|italic|W
      and Z bosons in association with two jets using the POWHEG method>,
      <hlink|<with|font-family|tt|1303.5447>|http://xxx.lanl.gov/abs/1303.5447>.

      <bibitem*|2><label|bib-Hamilton:2012np>K.<nbsp>Hamilton, P.<nbsp>Nason,
      and G.<nbsp>Zanderighi, <with|font-shape|italic|MINLO: Multi-Scale
      Improved NLO>, <with|font-shape|italic|JHEP>
      <with|font-series|bold|1210> (2012) 155,
      [<hlink|<with|font-family|tt|1206.3572>|http://xxx.lanl.gov/abs/1206.3572>].

      <bibitem*|3><label|bib-Nason:2013uba>P.<nbsp>Nason and C.<nbsp>Oleari,
      <with|font-shape|italic|Generation cuts and Born suppression in
      POWHEG>, <hlink|<with|font-family|tt|1303.3922>|http://xxx.lanl.gov/abs/1303.3922>.

      <bibitem*|4><label|bib-Melia:2011gk>T.<nbsp>Melia, P.<nbsp>Nason,
      R.<nbsp>Rontsch, and G.<nbsp>Zanderighi,
      <with|font-shape|italic|W<rsup|+>W<rsup|+> plus dijet production in the
      POWHEGBOX>, <with|font-shape|italic|Eur.Phys.J.>
      <with|font-series|bold|C71> (2011) 1670,
      [<hlink|<with|font-family|tt|1102.4846>|http://xxx.lanl.gov/abs/1102.4846>].
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|info-flag|detailed>
    <associate|par-hyphen|normal>
    <associate|preamble|false>
    <associate|sfactor|4>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|7|5>>
    <associate|auto-11|<tuple|6|?>>
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|3|1>>
    <associate|auto-4|<tuple|4|2>>
    <associate|auto-5|<tuple|4.1|2>>
    <associate|auto-6|<tuple|5|4>>
    <associate|auto-7|<tuple|6|4>>
    <associate|auto-8|<tuple|7|5>>
    <associate|auto-9|<tuple|7|5>>
    <associate|bib-Alioli:2010xd|<tuple|2|?>>
    <associate|bib-Alwall:2006yp|<tuple|1|5>>
    <associate|bib-Boos:2001cv|<tuple|2|6>>
    <associate|bib-Campbell:1999ah|<tuple|3|2>>
    <associate|bib-Campbell:2011bn|<tuple|4|2>>
    <associate|bib-Campbell:2013vha|<tuple|1|5>>
    <associate|bib-Dixon:1998py|<tuple|2|2>>
    <associate|bib-Frixione:2007nw|<tuple|1|?>>
    <associate|bib-Frixione:2007vw|<tuple|3|?>>
    <associate|bib-Frixione:2007zp|<tuple|3|6>>
    <associate|bib-Hamilton:2012np|<tuple|2|5>>
    <associate|bib-Mangano:1992jk|<tuple|4|6>>
    <associate|bib-Melia:2011gk|<tuple|4|5>>
    <associate|bib-Nason:1988xz|<tuple|5|6>>
    <associate|bib-Nason:1989zy|<tuple|6|6>>
    <associate|bib-Nason:2004rx|<tuple|7|6>>
    <associate|bib-Nason:2013uba|<tuple|3|5>>
    <associate|bib-noi|<tuple|1|2>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      Campbell:2013vha

      Hamilton:2012np

      Nason:2013uba

      Campbell:2013vha

      Melia:2011gk

      Melia:2011gk

      Campbell:2013vha
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Running
      the program> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Process
      specific input parameters> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Option
      that activate the new <with|font-family|<quote|tt>|POWHEG BOX>
      features> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|4.1<space|2spc>Parallel runs
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Files
      to inspect at the end of the run> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Shower
      and analysis> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|7<space|2spc>Examples>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>