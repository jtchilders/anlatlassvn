<TeXmacs|1.0.7.3>

<style|article>

<\body>
  <doc-data|<doc-title|The POWHEG-BOX-hvq manual>|>

  <section|Introduction>

  The POWHEG-hvq program can be used to generate heavy quark events (i.e.
  <math|t<wide|t|\<bar\>>>, <math|b<wide|b|\<bar\>>> and
  <math|c<wide|c|\<bar\>>> events) in hadronic collisions. It was first
  implemented in ref. <cite|Frixione:2007nw>, and then it was implemented
  again within the framework of the <with|font-family|tt|POWHEG BOX>
  <cite|Alioli:2010xd>. If you use this code, please quote
  <cite|Frixione:2007nw> and <cite|Alioli:2010xd>.

  This document describes the input parameters that are specific to this
  implementation. The parameters that are common to all
  <with|font-family|tt|POWHEG BOX> implementation are given in the
  <with|font-family|tt|manual-BOX.pdf> document, in the
  <with|font-family|tt|POWHEG-BOX/Docs> directory.

  <section|Generation of events>

  Do<next-line><with|font-family|tt|$ cd POWHEG-BOX/hvq><next-line><with|font-family|tt|$
  make pwhg_main><next-line>Then do (for example)<next-line><with|font-family|tt|$
  cd testrun-b-lhc><next-line><with|font-family|tt|$
  ../pwhg_main><next-line>At the end of the run, the file
  <with|font-family|tt|pwgevents.lhe> will contain 1000000 events for
  <math|b> pair production in the Les Houches format. In order to shower them
  with <with|font-family|tt|PYTHIA>:<next-line><with|font-family|tt|$ cd
  POWHEG-BOX/hvq><next-line><with|font-family|tt|$ make
  main-PYTHIA-lhef><next-line><with|font-family|tt|$ cd
  testrun-b-lhc><next-line><with|font-family|tt|$ ../main-PYTHIA-lhef>

  <section|Input parameters>

  Parameters in <with|font-family|tt|powheg.input> that are specific to heavy
  flavour production:<next-line><with|font-family|tt|qmass 4.75
  \ \ \ \ \ \ \ \ \ ! mass of heavy quark in GeV><next-line>

  The reference factorization and renormalization scales are taken by default
  equal to <math|<sqrt|p<rsub|T><rsup|2>+m<rsup|2><rsub|>>>, where
  <math|p<rsub|T>> is the transverse momentum of the heavy quark in the
  underlying Born configuration (i.e. before radiation).
  If<next-line><with|font-family|tt|fixedscale 1 \ \ \ \ \ \ \ \ ! use
  reference ren. and fact. scale = qmass><next-line>the reference scale is
  <math|m>

  <subsection|Top decays>

  \;

  If the quark being produced is a top, which is automatically determined by
  the program on the basis of the value of the mass, the top decay may be
  driven by the flag<next-line><next-line><with|font-family|tt|topdecaymode
  20000 ! an integer of 5 digits representing the decay
  mode>.<next-line><next-line>Top is assumed to go to a <math|b> and a
  <math|W>, with the <math|W> decaying according to a diagonal CKM matrix.The
  meaning of the token is the following: each digit represents the maximum
  number of the following particles \ in the (parton level) decay of the
  <math|t<wide|t|\<bar\>>> pair: <math|e<rsup|\<pm\>>>,
  <math|\<mu\><rsup|\<pm\>>>, <math|\<tau\><rsup|\<pm\>>>,
  <math|u<rsup|\<pm\>>>, <math|c<rsup|\<pm\>>>. Thus, for example, 20000
  means the <math|t\<rightarrow\>e<rsup|+>\<nu\><rsub|e> b>,
  <math|<wide|t|\<bar\>>\<rightarrow\>e<rsup|-><wide|\<nu\>|\<bar\>><rsub|e><wide|b|\<bar\>>>,
  22222 means all decays, 10011 means one goes into eletron or antielectron,
  and the other goes into any hadron, 00022 means fully hadronic, 00011 means
  fully hadronic with a single charm, 00012 fully hadronic with at least one
  charm. The value 0 means that the <math|t> and <math|<wide|t|\<bar\>>> are
  not decayed. Values that imply only one <math|t> decay (for example 10000)
  are not implemented consistently.

  If the flag <with|font-family|tt|semileptonic> is set to 1, only
  semileptonic decays are kept by the program.

  In case <with|font-family|tt|topdecaymode> is different from 0 more
  parameters are needed for the decay kinematics, and are used exclusively
  for decays<next-line><next-line><with|font-family|tt|tdec/wmass 80.4
  \ \ \ ! W mass for top decay<next-line>tdec/wwidth 2.141 \ ! W
  width<next-line>tdec/bmass 5 \ \ \ \ \ \ ! b quark mass in t
  decay<next-line>tdec/twidth 1.31 \ \ ! top width<next-line>tdec/sin2w 0.23
  \ \ \ ! Weinberg angle<next-line>tdec/elbranching 0.108 \ ! W electronic
  branching fraction<next-line>tdec/emass 0.00051 ! electron
  mass<next-line>tdec/mumass 0.1057 ! mu mass<next-line>tdec/taumass 1.777 !
  tau mass<next-line>tdec/dmass 0.100 \ \ ! d mass<next-line>tdec/umass 0.100
  \ \ ! u mass<next-line>tdec/smass 0.200 \ \ ! s mass<next-line>tdec/cmass
  1.5 \ \ \ \ ! charm mass<next-line>tdec/sin2cabibbo 0.051 \ ! sine of
  Cabibbo angle><next-line><next-line>If <with|font-family|tt|topdecaymode>
  is not set, or is set to 0, then top decay is not performed by
  <with|font-family|tt|POWHEG>. The event is passed to the Monte Carlo with
  undecayed tops, and the Shower program drives the decay. In this case, no
  spin correlations for the decay are included.

  <\bibliography|bib|JHEP|paper.bib>
    <\bib-list|1>
      <bibitem*|1><label|bib-Frixione:2007nw>S.<nbsp>Frixione, P.<nbsp>Nason,
      and G.<nbsp>Ridolfi, <with|font-shape|italic|A Positive-weight
      next-to-leading-order Monte Carlo for heavy flavour hadroproduction>,
      <with|font-shape|italic|JHEP> <with|font-series|bold|0709> (2007) 126,
      [<hlink|<with|font-family|tt|0707.3088>|http://xxx.lanl.gov/abs/0707.3088>].

      <bibitem*|2><label|bib-Alioli:2010xd>S.<nbsp>Alioli, P.<nbsp>Nason,
      C.<nbsp>Oleari, and E.<nbsp>Re, <with|font-shape|italic|A general
      framework for implementing NLO calculations in shower Monte Carlo
      programs: the POWHEG BOX>, <with|font-shape|italic|JHEP>
      <with|font-series|bold|1006> (2010) 043,
      [<hlink|<with|font-family|tt|1002.2581>|http://xxx.lanl.gov/abs/1002.2581>].
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|info-flag|detailed>
    <associate|preamble|false>
    <associate|sfactor|5>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|7|5>>
    <associate|auto-11|<tuple|6|?>>
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|3|1>>
    <associate|auto-4|<tuple|3.1|2>>
    <associate|auto-5|<tuple|3.1|2>>
    <associate|auto-6|<tuple|6|5>>
    <associate|auto-7|<tuple|7|5>>
    <associate|auto-8|<tuple|7|5>>
    <associate|auto-9|<tuple|7|5>>
    <associate|bib-Alioli:2010xd|<tuple|2|?>>
    <associate|bib-Alwall:2006yp|<tuple|1|5>>
    <associate|bib-Boos:2001cv|<tuple|2|6>>
    <associate|bib-Frixione:2007nw|<tuple|1|?>>
    <associate|bib-Frixione:2007vw|<tuple|3|?>>
    <associate|bib-Frixione:2007zp|<tuple|3|6>>
    <associate|bib-Mangano:1992jk|<tuple|4|6>>
    <associate|bib-Nason:1988xz|<tuple|5|6>>
    <associate|bib-Nason:1989zy|<tuple|6|6>>
    <associate|bib-Nason:2004rx|<tuple|7|6>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      Frixione:2007nw

      Alioli:2010xd

      Frixione:2007nw

      Alioli:2010xd
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Generation
      of events> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Input
      parameters> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|3.1<space|2spc>Top decays
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>