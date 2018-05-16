<TeXmacs|1.0.7.3>

<style|article>

<\body>
  <doc-data|<doc-title|Weighted events in the <with|font-family|tt|BOX>>>

  Two flags control the nature of the output in the
  <with|font-family|tt|BOX>: <with|font-family|tt|withnegweights> and
  <with|font-family|tt|ptsupp>. If neither of these flags is set, events are
  output with weight 1, i.e. the <with|font-family|tt|XWGTUP> variable is set
  to 1. The <with|font-family|tt|IDTUP> variable in the Les Houches interface
  is set to 3 in this case. The total cross section is stored in the
  <with|font-family|tt|xsecup> variable. Negative weighted events are
  neglected with this choice.

  If \ <with|font-family|tt|withnegweights> is set to 1 (true in our
  convention), negative weighted events are not discarded. The
  <with|font-family|tt|IDTUP> variable in the Les Houches interface is set to
  -4, and the weight of the event is set to its sign times the total cross
  section for positive weighted events plus the (absolute value) of the cross
  section for negative weighted events. In this way, the average value of
  <with|font-family|tt|XWGTUP> equals the real cross section, as required by
  the Les Houches convention when the <with|font-family|tt|IDTUP> variable is
  set to -4. The variable <with|font-family|tt|xsecup> always stores the real
  cross section.

  If the <with|font-family|tt|ptsupp> token is set, a suppression factor that
  depends upon the underlying Born configuration of each event is supplied
  with it. The cross section computed by the <with|font-family|tt|pwhg_main>
  program is in this case not valid. It is the integral of the cross section
  times the suppression factor. Events are generated using this ``fake''
  cross section, and thus are weighted with the inverse of the suppression
  factor. The <with|font-family|tt|IDTUP> variable in the Les Houches
  interface is set to -4. The weight of the event is in this case the sign,
  times the total cross section for positive weighted events plus the
  (absolute value) of the cross section for negative weighted events, times
  the inverse of the weighting factor. The weight factor is returned by the
  user routine \ <with|font-family|tt|born_suppression>, that can use the
  value of the <with|font-family|tt|ptsupp> token as a parameter to compute
  the suppression factor. Also in this case, the average value of the weight
  of the event is equal to the real cross section. This option can be active
  in conjunction with the <with|font-family|tt|withnegweights> flag.

  These flags have many uses. On one side, one might like to know where
  negative weighted events end up. The fact that they constitute a small
  fraction leaves us to worry that they may end up in some tiny tail of some
  important distribution. One may also prefer to work with negative weight in
  cases when getting rid of them requires high folding numbers (the
  <with|font-family|tt|foldcsi>, <with|font-family|tt|foldy> and
  <with|font-family|tt|foldphi> tokens), and thus has a high cost in computer
  time. The <with|font-family|tt|ptsupp> feature can be use to enhance a
  region of phase space (like a high <math|k<rsub|T>> tail) where it would be
  otherwise difficult to get high statistics. These features, however, become
  really useful for processes where the Born contribution itself is singular.
  The simplest examples are the <math|Z+jet> and the dijet production
  processes. Here we discuss <math|Z<rsub|>+jet>. The dijet case is fully
  analogous.

  \;

  \;

  <subsection|Generation cut and Born suppression factor>

  The <with|mode|math|Z+1j> process differs substantially from all processes
  previously implemented in <with|font-family|tt|POWHEG>, in the fact that
  the Born diagram itself is collinear and infrared divergent. In all
  previous implementations, the Born diagram was finite, and it was thus
  possible to generate an unweighted set of underlying Born configurations
  covering the whole phase space. In the present case, this is not possible,
  since they would all populate the very low transverse momentum region. Of
  course, this problem is also present in standard Shower Monte Carlo
  programs, where it is dealt with by generating the Born configuration with
  a cut <math|k<rsub|gen>> <math|k<rsub|gen>> on the transverse momentum of
  the <with|mode|math|Z> boson. After the shower, one must discard all events
  that fail some transverse momentum analysis cut <math|k<rsub|an>>
  <math|k<rsub|an>> in order to get a realistic sample. The analysis cut
  <math|k<rsub|an>> may be applied to the transverse momentum of the
  <with|mode|math|Z>, or to the hardest jet. We assume here, for sake of
  discussion that the analysis cut is applied to the <with|mode|math|Z>
  transverse momentum.

  Taking <with|mode|math|k<rsub|an>\<gtrsim\>k<rsub|gen>> is not enough to
  get a realistic sample. In fact, in an event generated at the Born level
  with a given <with|mode|math|k<rsub|T>\<less\>k<rsub|gen>>, the shower may
  increase the transverse momentum of the jet so that
  <with|mode|math|k<rsub|T><rsup|Z>\<gtr\>k<rsub|an>>. Thus, the generation
  cut, even if it is below the analysis cut, may reduce the number of events
  that pass the analysis cut. Of course, as we lower <math|k<rsub|gen>>
  keeping <math|k<rsub|an>> fixed, we will reach a point when very few events
  below <math|k<rsub|gen>> will pass the analysis cut <math|k<rsub|an>>. In
  fact, generation of radiation with transverse momentum larger than
  <math|k<rsub|gen>> is strongly suppressed in <with|font-family|tt|POWHEG>,
  and, in turn, radiation from subsequent shower is required to be not harder
  than the hardest radiation of <with|font-family|tt|POWHEG>. Thus, given the
  fact that we want to generate a sample with a given <math|k<rsub|an>> cut,
  we should choose <math|k<rsub|gen>> small enough, so that the final sample
  remains substantially the same if <math|k<rsub|gen>> is lowered even
  further.

  A second option for the implementation of processes with a divergent Born
  contribution is also available. It requires that we generate weighed
  events, rather than unweighted ones. This is done by using a suppressed
  cross section for the generation of the underlying Born configurations:

  <\equation>
    <wide|B|\<bar\>><rsub|<with|mode|text|font-family|rm|supp>>=<wide|B|\<bar\>>\<times\>F(k<rsub|T>),
  </equation>

  wher <with|mode|math|<wide|B|\<bar\>>> is the inclusive NLO cross section
  at fixed underlying Born variables, and <with|mode|math|k<rsub|T>> is the
  transverse momentum of the vector boson in the underlying Born
  configuration. In this way <with|mode|math|<wide|B|\<bar\>><rsub|<with|mode|text|font-family|rm|supp>>>
  is integrable, and one can use it to generate underlying Born
  configurations according to its value. The generated event, however, should
  be given a weight <with|mode|math|1/F(k<rsub|T>)> rather than a constant
  one, in order to compensate for the initial <with|mode|math|F(k<rsub|T>)>
  suppression factor. With this method, events do not concentrate in the low
  <with|mode|math|k<rsub|T>> region. However, their weight in the low
  <with|mode|math|k<rsub|T>> region becomes divergent. After shower, if one
  imposes the analysis cut, one gets a finite cross section, since it is
  unlikely that events with small transverse momentum at the Born level may
  pass the analysis cut after shower. In fact, shower transverse momenta
  larger than the one present in the initial Born process must be suppressed
  in the Monte Carlo generator.

  In recent <with|font-family|tt|POWHEG BOX> revisions, both methods can be
  implemented at the same time. We wanted i fact to be able to implement the
  following three options:

  <\itemize>
    <item>Generate events using a transverse momentum generation cut.

    <item>Generate events using a Born suppression factor, and a small
    transverse momentum cut, just enough to avoid unphysical values of the
    strong coupling constant and of the factorization scale that appears in
    the parton density functions.

    <item>Apply a Born suppression factor, and set the transverse momentum
    cut to zero. In this case the program cannot be used to generate events.
    It can be used, however, to produce NLO fixed order distributions,
    provided the renormalization and factorization scales are set in such a
    way that they remain large enough even at small
    <with|mode|math|k<rsub|T><rsup|Z>>. This feature is only used for the
    generation of fixed order distributions.
  </itemize>

  The generation cut is activated by setting the token
  <with|font-family|tt|bornktmin> to the desired value in the
  <with|font-family|tt|powheg.input> file. The Born suppression is activated
  by setting the token <with|font-family|tt|ptsupp> to a positive real value.
  The process-specific subroutine <with|font-family|tt|born_suppression> sets
  the suppression factor to <with|mode|math|k<rsub|T><rsup|2>/(k<rsub|T><rsup|2>+<with|math-font-family|mt|ptsupp><rsup|2>)>.
  If <with|font-family|tt|psupp> is negative, the suppression factor is set
  to 1.

  The need of a transverse momentum cut is not only a technical issue. The
  NLO calculation of <with|mode|math|Z+1j> production holds only if the
  transverse momentum of the <with|mode|math|Z> is not too small. In fact, as
  the <math|k<rsub|T>> decreases, large Sudakov logarithms arise in the NLO
  correction, and the value of the running coupling increases, up to the
  point when the cross section at fixed order becomes totally unreliable.
  These large logarithms should all be resummed in order to get a sensible
  answer in this region. In the <with|font-family|tt|POWHEG> implementation
  of <with|mode|math|Z> production, in fact, these logarithms are all
  resummed. It is clear then that some sort of merging between the
  <with|mode|math|Z+1j> and the <with|mode|math|Z> production processes
  should be performed at relatively small transverse momentum, in order to
  properly deal with these large logarithms. In the present work we will not
  attempt to perform such merging, that we leave for future publications. We
  will simply remember, when looking at our results, that we expect to get
  unphysical distributions when the <with|mode|math|Z> transverse momentum is
  too small, and we will discuss this fact in a more quantitative way.

  In the <with|font-family|tt|POWHEG> approach, negative weighted events can
  only arise if one is approaching a region where the NLO computation is no
  longer feasible. In our studies for the <with|mode|math|Z+1j> process we
  approach this region at small transverse momentum. In order to better see
  what happens there, rather than neglecting negative weights (that is the
  default behaviour of the <with|font-family|tt|POWHEG BOX>), we have
  introduced a new feature in the program, that allows one to track also the
  negative weighted events. This feature is activated by setting the token
  <with|font-family|tt|withnegweights> to 1 (true). If
  \ <with|font-family|tt|withnegweights> is set to 1, events with negative
  weight can thus appear in the Les Houches event file. While we normally set
  the <with|font-family|tt|IDWTUP> flag in the Les Houches interface to 3, in
  this case we set it to -4. With this flag, the SMC is supposed to simply
  process the event, without taking any other action. Furthermore, the
  <with|font-family|tt|XWGTUP> (Les Houches) common block variable is set by
  the <with|font-family|tt|POWHEG-BOX> to the sign of the event times the
  integral of the absolute value of the cross section, in such a way that its
  average equals the true total cross section.

  Notice that, if <with|font-family|tt|withnegweights> is set and a Born
  suppression factor is also present, the events will have variable
  <with|font-family|tt|XWGTUP> of either signs. In this case
  <with|font-family|tt|XWGTUP> is set to the sign of the event, times the
  absolute value of the cross section, divided by the suppression factor
  <with|font-family|tt|ptsupp>. Also in this case the average value of
  <with|font-family|tt|XWGTUP> coincides withe the true total cross section.
  We preferred not to use the option -3 in case of signed events with
  constant absolute value. This option is advocated by the Les Houches
  interface precisely in such cases. However, the Les Houches interface does
  not provide a standard way to store the integral of the absolute value of
  the cross section, that would be needed to compute correctly the weight of
  the event. In fact, the <with|font-family|tt|XSECUP> variable is reserved
  for the true total cross section. More specifically, if we have
  <with|mode|math|N> events of either sign, they should be weighted with the
  sum of the positive plus the absolute value of the negative part of the
  cross section, in such a way that

  <\equation>
    <big|sum><rsub|i=1><rsup|N>W<rsub|i><left|(>\<sigma\><rsub|(+)>+\|\<sigma\><rsub|(-)>\|<right|)>=N<left|(>\<sigma\><rsub|(+)>-\|\<sigma\><rsub|(-)>\|<right|)>=<space|0.75spc>N\<sigma\><rsub|<with|mode|text|font-family|rm|NLO>>,
  </equation>

  (where <with|mode|math|W<rsub|i>> are the sign of the event
  <with|mode|math|\<pm\>1>), because

  <\equation>
    <frac|<big|sum><rsub|i>W<rsub|i>|N>=<frac|<left|(>\<sigma\><rsub|(+)>-\|\<sigma\><rsub|(-)>\|<right|)>|<left|(>\<sigma\><rsub|(+)>+\|\<sigma\><rsub|(-)>\|<right|)>><space|0.75spc>.
  </equation>

  Weighted events are also useful if one wants to generate a homogeneous
  sample from relatively low up to very high transverse momenta. It is
  convenient in this case to pick a very large <with|font-family|tt|ptsupp>
  value, of the order of the maximum transverse momentum one is interested
  in. The large momentum region will be more populated in this way.
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      Boos:2001cv

      Alwall:2006yp
    </associate>
    <\associate|toc>
      <with|par-left|<quote|1.5fn>|1<space|2spc>Generation cut and Born
      suppression factor <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>