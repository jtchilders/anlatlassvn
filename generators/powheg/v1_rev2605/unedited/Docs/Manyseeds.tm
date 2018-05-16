<TeXmacs|1.0.7.3>

<style|generic>

<\body>
  <doc-data|<doc-title|The <with|font-family|tt|manyseeds> flag>>

  If the <with|font-family|tt|manyseeds> flag is set to 1, the
  <with|font-family|tt|BOX> looks for a file named
  <with|font-family|tt|<with|font-family|tt|pwgseeds>.dat>. The program stops
  if this file is not found. If the file is found the program asks for an
  integer <with|font-family|tt|j> for input. Let us say that this integer is
  17, for sake of clarity. This integer indicates the line in
  <with|font-family|tt|pwgseed.dat> where a random number seed is stored.
  This feature is used for running the <with|font-family|tt|BOX> on several
  nodes of a cluster. In order to do this, follow the steps

  <\enumerate-numeric>
    <item>Prepare a <with|font-family|tt|powheg.input> file with the
    <with|font-family|tt|manyseeds> flag is set to 1. Set the number of
    events <with|font-family|tt|nev> to 0.

    <item>Prepare a file <with|font-family|tt|pwgseeds.dat>, containing a
    sequence (one per line) of different random number seeds. For example:
    first line 1, second line 2, etc. (but any number will do).

    <item>Run the <with|font-family|tt|pwhg_main> program. It will ask an
    integer for input. Input an integer. That integer is the line number of
    the random seed to be used for the current run. Assuming that the number
    17 is given as input to the <with|font-family|tt|pwhg_main> program, the
    run will produce files named <with|font-family|tt|pwgxgrid.dat> and
    <with|font-family|tt|pwggrid-0017.dat>. The
    <with|font-family|tt|pwhg_main> program can be run with different
    integers as input. Each run can be sent, for example, to a different node
    of a cluster.

    At the end of this step, a bunch of <with|font-family|tt|pwggrid-[????].dat>,
    <with|font-family|tt|pwgubound-[????].dat> files ans
    <with|font-family|tt|pwgNLO-[????].top> files will be present in the run
    directory. The <with|font-family|tt|pwgNLO-[????].top> files are
    statistically independent topdrawer histograms. They can be combined to
    provide a higher statistics NLO analysis of the current analysis
    routines.

    If the subsequent runs are sent after the file
    <with|font-family|tt|pwgxgrid.dat> was already produced, and if the flag
    <with|font-family|tt|><with|font-family|tt|use-old-grid> is set to one,
    the importance sampling grid will be loaded from the
    \ <with|font-family|tt|pwgxgrid.dat>. Otherwise it will be recreated.

    <item>Now set the <with|font-family|tt|nev> token to a given number, make
    sure that the flag <with|font-family|tt|><with|font-family|tt|use-old-grid>
    is set to one, and run a bunch of copies of the
    <with|font-family|tt|pwhg_main> program, each with the same integers as
    input. The program will now load all the
    <with|font-family|tt|pwggrid-[????].dat> and
    <with|font-family|tt|pwgubound-[????].dat> that it can find, and combines
    them adequately, assuming that they are all statistically independent,
    and will start to generate events. The events will be in files
    <with|font-family|tt|pwgevent-[????].lhe>, and they will all be
    independent statistically.
  </enumerate-numeric>

  The sequence above is a quite simple two step procedure. It is useful,
  however, to better clarify the logic that the <with|font-family|tt|BOX>
  follows in this procedure.

  There are 3 steps in the initialization phase of
  <with|font-family|tt|POWEG>. First of all, an importance sampling grid is
  determined. Let us call this stage ISG (Importance Sampling Grid). The
  second step is the calculation of the integrals, and the determination of
  an upper bounding envelope for the <math|<wide|B|~>> function, to be used
  for the generation of underlying Born configurations. We call this stage
  the UBB (Uppe Bounds for underlying Born). As a third step, the upper bound
  normalization for radiation is \ determined. We call this stage UBR (Upper
  Bounds for Radiation).

  \ First of all, we remark that, if the \ <with|font-family|tt|><with|font-family|tt|use-old-grid>
  flag is not set to 1, no grid file is loaded. Similarly, if the
  <with|font-family|tt|><with|font-family|tt|use-old-ubound> flag is not set
  to 1, no ubound file is loaded. In other words, if these flags are set,
  looking for a corresponding file will always yield a negative result. The
  reader should keep this in mind when reading the following procedure. If
  the flag <with|font-family|tt|manyseeds> is set, the
  <with|font-family|tt|pwhg_main> program asks for an integer. We will call
  <with|font-family|tt|cj> this integer, between 0 and 9999. We denote with
  <with|font-shape|right|[cj]> the corresponding string of four digits
  (leading digits are set to 0; thus if <with|font-family|tt|ic=1>
  <with|font-family|tt|[cj]=0001>. We will denote as
  <with|font-family|tt|[????]> any four digit string. The logic of grid
  loading in the BOX is as follows:

  <\enumerate-numeric>
    <item>If a file <with|font-family|tt|pwggrid.dat> exist and is
    consistent, this file is loaded, and steps ISG and UBB are skipped (go to
    <reference|step:UBR>).

    <item>If the above fails, if the <with|font-family|tt|manyseeds> flag is
    not set, or if it is set and a file named
    <with|font-family|tt|pwggrid-[cj].dat> already exists and is consistent,
    all files of the form <with|font-family|tt|pwggrid-[????].dat> are loaded
    and suitably combined, and the steps ISG and UBB are skipped (go to
    <reference|step:UBR>).

    <item>If the above fails, if a file <with|font-family|tt|pwgxgrid.dat>
    \ exist and is consistent, this file is loaded, and steps ISG is skipped
    (go to <reference|step:UBB><active*|>).

    <item><label|step:ISG>Step ISG is performed. The resulting grid is stored
    in the file <with|font-family|tt|pwgxgrid.dat>. This step, whether the
    <with|font-family|tt|manyseeds> flag is set or not, is performed using
    the default initial seed value (i.e. not the seed found at the
    <with|font-family|tt|cj> line of the <with|font-family|tt|pwgseeds.dat>
    file). In this way, all copies of the program being run will use the same
    importance sampling greed. This is mandatory if we want to combine
    results.

    <item><label|step:UBB>Step UBB is performed. If the
    <with|font-family|tt|manyseeds> flag is set, this step is performed using
    the seed found at the <with|font-family|tt|cj> line of the
    <with|font-family|tt|pwgseeds.dat> file, and the result is storedin the
    file <with|font-family|tt|pwggrid-[cj].dat>. Otherwise, the current seed
    value is used, and the result is stored in a file named
    <with|font-family|tt|pwggrid.dat>.

    <item><label|step:UBR>If a file named <with|font-family|tt|pwgubound.dat>
    exists and is consistent it is loaded. The UBR step is skipped (goto
    <reference|step:GENEV>).

    <item>If a file named <with|font-family|tt|pwgubound.dat> exists and is
    consistent it is loaded. The UBR step is skipped (goto
    <reference|step:GENEV>).

    <item>If the above fails, if the <with|font-family|tt|manyseeds> flag is
    not set, or if it is set and a file named
    <with|font-family|tt|pwgubound-[cJ].dat> exists and is consistent, all
    files with names of the form <with|font-family|tt|pwgubound-[????].dat>
    are loaded and combined. The UBR step is skpped (go to
    <reference|step:GENEV>).

    <item>The UBR step is performed. \ If the <with|font-family|tt|manyseeds>
    flag is not set, the result is stored in a file named
    <with|font-family|tt|pwgubound.dat>. Otherwise, it is stored in the file
    <with|font-family|tt|pwgubound-[cj].dat>.

    <item><label|step:GENEV>Now <with|font-family|tt|nev> events are
    generated. If the <with|font-family|tt|manyseeds> flag is not set, the
    result is stored in a file named <with|font-family|tt|pwgevents.lhe>.
    Otherwise, it is stored in the file <with|font-family|tt|pwgevents-[cj].dat>.
  </enumerate-numeric>

  This logic has the purpose to allow several possible combinations of
  actions. For example, one can use grids generated in parallel runs to
  produce events without using the <with|font-family|tt|manyseeds> flag. Or
  one can use grids generated without the <with|font-family|tt|manyseeds>
  flag for generating events in parallel with the
  <with|font-family|tt|manyseeds> flag set.
</body>

<\initial>
  <\collection>
    <associate|sfactor|4>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|step:GENEV|<tuple|10|?>>
    <associate|step:ISG|<tuple|4|?>>
    <associate|step:UBB|<tuple|5|?>>
    <associate|step:UBR|<tuple|6|?>>
  </collection>
</references>