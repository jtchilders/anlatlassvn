<TeXmacs|1.0.7.3>

<style|generic>

<\body>
  Write the <math|n+1> body phase space in factorized form

  <\equation>
    d\<Phi\><rsub|n+1>=<frac|d M<rsub|rec><rsup|2>|2\<pi\>>
    <frac|d<rsup|3>k<rsub|n+1>|2k<rsub|n+1><rsup|0>(2\<pi\>)<rsup|3>>
    <frac|d<rsup|3>k<rsub|n>|2k<rsub|n><rsup|0>(2\<pi\>)<rsup|3>>
    <frac|d<rsup|3>k<rsub|rec>|2k<rsub|rec><rsup|0>(2\<pi\>)<rsup|3>>(2\<pi\>)<rsup|4>\<delta\><rsup|4>(q-k<rsub|n+1>-k<rsub|n>-k<rsub|rec>)
    \ d\<Phi\><rsub|rec> ,
  </equation>

  where <math|d\<Phi\><rsub|rec>> is the phase space of the system recoiling
  agains <math|k<rsub|n+1>+k<rsub|n>>. The three body phase space part can be
  written in the Daliz form:

  <\equation*>
    <frac|d<rsup|3>k<rsub|n+1>|2k<rsub|n+1><rsup|0>(2\<pi\>)<rsup|3>>
    <frac|d<rsup|3>k<rsub|n>|2k<rsub|n><rsup|0>(2\<pi\>)<rsup|3>>
    <frac|d<rsup|3>k<rsub|rec>|2k<rsub|rec><rsup|0>(2\<pi\>)<rsup|3>>(2\<pi\>)<rsup|4>\<delta\><rsup|4>(q-k<rsub|n+1>-k<rsub|n>-k<rsub|rec>)=
  </equation*>

  <\equation*>
    <frac|d<rsup|3>k<rsub|n+1>|2k<rsub|n+1><rsup|0>(2\<pi\>)<rsup|3>>
    <frac|d<rsup|3>k<rsub|n>|2k<rsub|n><rsup|0>(2\<pi\>)<rsup|3>>
    2\<pi\>\<delta\>((q-k<rsub|n+1>-k<rsub|n>)<rsup|2>-M<rsub|rec><rsup|2>)=
  </equation*>

  <\equation*>
    <frac|d\<Omega\><rsub|n+1>|4(2\<pi\>)<rsup|6>> \ k<rsub|n+1> d
    k<rsub|n+1><rsup|0> \ \ k<rsub|n> d k<rsub|n><rsup|0> \ d cos\<theta\>
    d\<phi\> \ 2\<pi\> \<delta\>((q-k<rsub|n+1>-k<rsub|n>)<rsup|2>-M<rsub|rec><rsup|2>)
    .
  </equation*>

  Now

  <\equation*>
    (q-k<rsub|n+1>-k<rsub|n>)<rsup|2>-M<rsub|rec><rsup|2>=q<rsup|2>+m<rsub|n><rsup|2>-M<rsub|rec><rsup|2>-2q<rsup|0>(k<rsub|n+1><rsup|0>+k<rsub|n><rsup|0>)+2k<rsub|n+1><rsup|0>k<rsub|n><rsup|0>+2cos\<theta\>
    k<rsub|n+1>k<rsub|n>,
  </equation*>

  so the <math|d cos\<theta\>> integration yields a
  <math|1/(2k<rsub|n+1>k<rsub|n>)> factor, yielding

  <\equation>
    <frac|d\<Omega\><rsub|n+1>|8(2\<pi\>)<rsup|5>> \ \ d k<rsub|n+1><rsup|0>
    \ \ \ d k<rsub|n><rsup|0> \ \ \ \ d\<phi\> .
  </equation>

  The orientation can be taken relative to any of the three bodies. So, it
  might as well be written

  <\equation>
    <frac|d\<Omega\>|8(2\<pi\>)<rsup|5>> \ \ d k<rsub|n+1><rsup|0> \ \ \ d
    k<rsub|n><rsup|0> \ \ \ \ d\<phi\> ,
  </equation>

  where <math|\<Omega\>> is the direction of <math|k<rsub|rec>>, and
  <math|\<phi\>> is the azimuth of <math|k<rsub|n>> or <math|k<rsub|n+1>>
  relative to <math|k<rsub|rec>>.

  The underlying Born phase space can be factorized into a two body phase
  space times the phase space of the recoil system

  <\equation>
    d<wide|\<Phi\>|\<bar\>><rsub|n>=<frac|d M<rsub|rec><rsup|2>|2\<pi\>>
    \ <frac|d<rsup|3><wide|k|\<bar\>><rsub|n>|2<wide|k|\<bar\>><rsub|n><rsup|0>(2\<pi\>)<rsup|3>>
    <frac|d<rsup|3>k<rsub|rec>|2k<rsub|rec><rsup|0>(2\<pi\>)<rsup|3>>
    (2\<pi\>)<rsup|4>\<delta\><rsup|4>(q-<wide|k|\<bar\>><rsub|n>-k<rsub|rec>)
    d\<Phi\><rsub|rec>.
  </equation>

  We have

  <\equation>
    <frac|d<rsup|3><wide|k|\<bar\>><rsub|n>|2<wide|k|\<bar\>><rsub|n><rsup|0>(2\<pi\>)<rsup|3>>
    <frac|d<rsup|3>k<rsub|rec>|2k<rsub|rec><rsup|0>(2\<pi\>)<rsup|3>>(2\<pi\>)<rsup|4>\<delta\><rsup|4>(q-<wide|k|\<bar\>><rsub|n>-k<rsub|rec>)
    =<frac|d<rsup|3><wide|k|\<bar\>><rsub|n>|2<wide|k|\<bar\>><rsub|n><rsup|0>(2\<pi\>)<rsup|3>>
    2\<pi\>\<delta\><left|(>(q-<wide|k|\<bar\>><rsub|n>)<rsup|2>-M<rsub|rec><rsup|2><right|)>=<frac|d\<Omega\>|32\<pi\><rsup|2>>
    <frac|2<wide|k|\<bar\>><rsub|n>|q<rsup|0>> .
  </equation>

  We wish to express the phase space in terms of the underlying Born phase
  space, and the radiation variables, that we take equal to
  <math|k<rsup|0><rsub|n+1>>, <math|k<rsup|0><rsub|n>> and <math|\<phi\>>. In
  other words, we must identify

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|d\<Omega\>|8(2\<pi\>)<rsup|5>> \ \ d
    k<rsub|n+1><rsup|0> \ \ \ d k<rsub|n><rsup|0> \ \ \ \ d\<phi\> \ <frac|d
    M<rsub|rec><rsup|2>|2\<pi\>> \ \ d\<Phi\><rsub|rec>>|<cell|=>|<cell|J d
    k<rsub|n+1><rsup|0> \ \ \ d k<rsub|n><rsup|0> \ d\<phi\>
    d<wide|\<Phi\>|\<bar\>><rsub|n>>>|<row|<cell|>|<cell|=>|<cell|J d
    k<rsub|n+1><rsup|0> \ \ \ d k<rsub|n><rsup|0> d\<phi\>
    <frac|d\<Omega\>|32\<pi\><rsup|2>> <frac|2<wide|k|\<bar\>><rsub|n>|q<rsup|0>>
    \ <frac|d M<rsub|rec><rsup|2>|2\<pi\>>
    \ \ d\<Phi\><rsub|rec><eq-number>>>>>
  </eqnarray*>

  Cancelling common factors, we get

  <\equation>
    \ \ \ J=<frac|1|(2\<pi\>)<rsup|3>> <frac|q<rsup|0>|2<wide|k|\<bar\>><rsub|n>>
    ,
  </equation>

  and

  <\equation>
    d\<Phi\><rsub|n+1>=J d k<rsub|n+1><rsup|0> \ \ \ d k<rsub|n><rsup|0>
    d\<phi\> \ d<wide|\<Phi\>|\<bar\>><rsub|n> .
  </equation>

  <section|Mapping to <math|><with|mode|math|k<rsub|n+1><rsup|0>> to
  <with|mode|math|y>>

  It is more convenient to introduce the following coordinates. It turns out
  that <math|k<rsup|0><rsub|n+1>>, <math|k<rsub|n><rsup|0>> and
  <math|k<rsub|rec><rsup|0>> live in a convex Dalitz
  domain.<\float|float|tbh>
    <big-figure|<postscript|dalitz.eps|*5/8|||||>|<label|fig:dalitz>Dalitz
    region for <math|M<rsub|rec>= 0.3>, <math|m=0.2> .>
  </float> The boundary of that domain is defined by the configuration
  regions where the three vectors <math|<wide|k|\<vect\>><rsub|n+1>>,
  <math|<wide|k|\<vect\>><rsub|n>> and <math|<wide|k|\<vect\>><rsub|rec>> lie
  on the same line. The point at <math|k<rsub|n+1>=0> certainly belongs to
  that domain. We easily find that when <math|k<rsub|n+1>=0> we have

  <\equation>
    k<rsub|n><rsup|0>=<wide|k|^><rsub|n><rsup|0>\<equiv\><frac|q<rsup|2>+m<rsup|2>-M<rsub|rec><rsup|2>|2q>
    , \ \ \ \ \ \ \ \ \ \ \ k<rsup|0<rsub|>><rsub|rec>=<wide|k|^><rsup|0><rsub|rec>\<equiv\><frac|q<rsup|2>-m<rsup|2>+M<rsub|rec><rsup|2>|2q>
    .
  </equation>

  Since the Dalitz domain is convex, it should be possible to parametrize it
  as a function of two parameters, <math|k<rsup|0><rsub|n+1>> itself and
  <math|z>, with

  <\equation>
    k<rsup|0><rsub|n>=<wide|k|^><rsup|0><rsub|n>-z k<rsub|n+1>
    .<label|eq:kfromy>
  </equation>

  This corresponds, in Fig. <reference|fig:dalitz>, to parametrize the point
  <math|X> by giving the \ <math|k<rsup|0><rsub|n+1>> value, and the tangent
  <math|y> of the angle <math|<wide|X O A|^>>. Thus, the point
  <math|k<rsub|n><rsup|0>=><with|mode|math|<wide|k|^><rsup|0><rsub|n>>,
  <math|><with|mode|math|k<rsub|n+1>=0> belongs to the Dalitz region for all
  values of <math|z>. For a given value of <math|z>, there is then a maximum
  value of <math|k<rsub|n+1>>, such that the point is on the boundary of the
  Dalitz region (the point <math|Y> in the figure). It is characterized by
  the condition

  <\equation>
    \|<wide|k|\<vect\>><rsub|n+1>\|\<pm\>\|<wide|k|\<vect\>><rsub|n>\|\<pm\>\|<wide|k|\<vect\>><rsub|rec>\|=0
    ,<label|eq:triangleeq>
  </equation>

  that has to hold for at least one sign combination. Eq.
  (<reference|eq:triangleeq>) corresponds to the boundary of a triangular
  inequality. It is solved by squaring

  <\equation>
    (\|<wide|k|\<vect\>><rsub|n+1>\|\<pm\>\|<wide|k|\<vect\>><rsub|n>\|)<rsup|2>=k<rsub|n+1><rsup|2>+<wide|k|\<vect\>><rsub|n><rsup|2>
    \<pm\>2k<rsub|n+1>\|<wide|k|\<vect\>><rsub|n>\|=<wide|k|\<vect\>><rsub|rec><rsup|2>
    ,
  </equation>

  from which it follows again

  <\equation>
    <left|(>k<rsup|2><rsub|n+1>+<wide|k|\<vect\>><rsup|2><rsub|n>-<wide|k|\<vect\>><rsup|2><rsub|rec><right|)><rsup|2>=4k<rsub|n+1><rsup|2><wide|k|\<vect\>><rsub|n><rsup|2>
    .
  </equation>

  We can now use <math|<wide|k|\<vect\>><rsup|2><rsub|n>=k<rsup|0><rsub|n><rsup|2>-m<rsup|2>>,
  <math|<wide|k|\<vect\>><rsub|rec><rsup|2>=(q-k<rsup|0><rsub|n>-k<rsup|0><rsub|n+1>)<rsup|2>-M<rsub|rec><rsup|2>>,
  and eq. (<reference|eq:kfromy>), and we obtain the equation for
  <math|k<rsub|n+1>> of the form

  <\equation>
    4k<rsub|n+1><rsup|2><left|(>2k<rsub|n+1> q \ z(1-z)+q<rsup|2>z<rsup|2>-2
    q <wide|k|^><rsup|0><rsub|rec> z+M<rsub|rec><rsup|2><right|)>=0 .
  </equation>

  This yields a double solution <math|k<rsub|n+1>=0>, that we already knew,
  and

  <\equation>
    k<rsub|n+1>=-<frac|q<rsup|2>z<rsup|2> -2 q <wide|k|^><rsup|0><rsub|rec>
    z+M<rsub|rec><rsup|2>|2q z(1-z)> ,<label|eq:knp1max>
  </equation>

  which is the seeked maximum value of <math|k<rsub|n+1>>. The numerator on
  the right hand side of eq. (<reference|eq:knp1max>) vanishes for

  <\equation>
    z<rsub|1/2>=<left|(><wide|k|^><rsup|0><rsub|rec>\<pm\><sqrt|(<wide|k|^><rsup|0><rsub|rec>)<rsup|2>-M<rsub|rec><rsup|2>><right|)>/q
    .
  </equation>

  These correspond to the maximum and minimum <math|z> values allowed. The
  lines at fixed <math|z=z<rsub|1/2>> are the lower/upper tangent to the
  Dalitz region from the point <math|O>, shown in the figure.

  We now define

  <\equation>
    k<rsub|n+1>=<frac|\<xi\>q|2>, \ \ \ \ \ \ \ \ z=z<rsub|2>-(z<rsub|2>-z<rsub|1>)(y+1)/2
    ,
  </equation>

  where now <math|z> and <math|y> will play the role of the FKS variables.
  Notice that <math|y=1> corresponds to <math|z=z<rsub|1>>, which is the
  lower tangent in the plot. The upper tangent corresponds to
  <math|k<rsub|n>> nearly constant, i.e. to <math|k<rsub|n>> recoiling
  against <math|k<rsub|n+1>> and <math|k<rsub|rec>>, that are parallel. The
  lower tangent corresponds instead to <math|k<rsub|rec>> recoiling against
  <math|k<rsub|n+1>, k<rsub|n>>, that are collinear. Thus, <math|y=1>
  corresponds to the mass singularity, which is what we want.

  <section|Transverse momentum definition>

  We need to find an appropriate definition of the transverse momentum in
  this case. The commonly adopted one in the massless case is not appropriate
  now, because, in the soft limit, the maximum virtuality available for the
  emitted gluon is no longer limited by the real transverse momentum. The
  singularity in the propagator of the massive quark has the structure

  <\equation>
    <frac|1|(k+p)<rsup|2>-m<rsup|2>> ,
  </equation>

  where we have called now for simplicity <math|k<rsub|n+1>\<rightarrow\>k>,
  <math|k<rsub|n>\<rightarrow\>p>. We must figure out the maximum virtuality
  that the emitted gluon can achieved in order not to perturb the
  corresponding singularity. We can rewrite

  <\equation>
    (k+p)<rsup|2>-m<rsup|2>=2k\<cdot\>p+k<rsup|2>=2(p<rsup|0>k<rsup|0>-<wide*|p|\<wide-bar\>>k<rsup|p>)+k<rsup|2>
    ,
  </equation>

  where <math|k<rsup|p>> is the component of <math|k> parallel to <math|p>.
  We can now rewrite

  <\eqnarray*>
    <tformat|<table|<row|<cell|2(p<rsup|0>k<rsup|0>-<wide*|p|\<wide-bar\>>
    k<rsup|p>)>|<cell|=>|<cell|(p<rsup|0>+<wide*|p|\<wide-bar\>>)(k<rsup|0>-k<rsup|p>)+(p<rsup|0>-<wide*|p|\<wide-bar\>>)(k<rsup|0>+k<rsup|p>)>>|<row|<cell|>|<cell|=>|<cell|(p<rsup|0>+<wide*|p|\<wide-bar\>>)<frac|k<rsub|\<bot\>><rsup|2>+k<rsup|2>|k<rsup|0>+k<rsup|p>>+<frac|m<rsup|2>|p<rsup|0>+<wide*|p|\<wide-bar\>>>(k<rsup|0>+k<rsup|p>)
    .<eq-number>>>>>
  </eqnarray*>

  We see that, as in the massless limit case,
  <math|k<rsup|2>\<less\>k<rsup|2><rsub|\<perp\>>>, but now also

  <\equation>
    (p<rsup|0>+<wide*|p|\<wide-bar\>>)<frac|k<rsup|2>|k<rsup|0>+k<rsup|p>>\<less\><frac|m<rsup|2>|p<rsup|0>+<wide*|p|\<wide-bar\>>>(k<rsup|0>+k<rsup|p>)
    \ \ \<Longrightarrow\>k<rsup|2>\<less\><frac|m<rsup|2>|(p<rsup|0>+<wide*|p|\<wide-bar\>>)<rsup|2>>(k<rsup|0>+k<rsup|p>)<rsup|2>.
  </equation>

  More simply, we can write

  <\equation>
    (p<rsup|0>+<wide*|p|\<wide-bar\>>)<frac|k<rsup|2>|k<rsup|0>+k<rsup|p>>\<less\>2k\<cdot\>p,
    \ \ \ \ \ \ or \ \ \ \ \ \ \ k<rsup|2>\<less\><frac|k<rsup|0>+k<rsup|p>|p<rsup|0>+<wide*|p|\<wide-bar\>>>
    2k\<cdot\>p ,
  </equation>

  where the expression on the right hand side coincides with the transverse
  momentum in the massless limit. In fact, in the massless limit it becomes

  <\equation>
    <frac|k<rsup|0>(1+cos \<theta\>)|2p<rsup|0>>2k\<cdot\>p \ \ .
  </equation>

  In the transverse momentum definition used in the
  \ <with|font-family|tt|POWHEG> <with|font-family|tt|BOX> the factor
  <math|(1+cos \<theta\>)/2> is usually dropped. Thus, even here we can use a
  transverse momentum definition

  <\equation>
    K<rsub|T><rsup|2>=<frac|2k<rsup|0>|(p<rsup|0>+<wide*|p|\<wide-bar\>>)>
    2k\<cdot\>p \ \ ,
  </equation>

  that in the massless limit becomes

  <\equation*>
    <frac|k<rsup|0>|p<rsup|0>> 2 k\<cdot\>p=2(k<rsup|0>)<rsup|2>(1-cos
    \<theta\>)
  </equation*>

  This corresponds exactly to the <with|font-family|tt|POWHEG>
  <with|font-family|tt|BOX> definition. We can further simplify our
  expression by taking <math|<wide*|p|\<wide-bar\>>\<rightarrow\>p<rsup|0>>,
  and thus adopt the formula

  <\equation>
    K<rsub|T><rsup|2>=<frac|k<rsup|0>|p<rsup|0>> 2 k\<cdot\>p \ \ .
  </equation>

  With the Dalitz variables defined above, both <math|k\<cdot\>p> and
  <math|p<rsup|0>> have a simple expression. We have

  <\equation>
    p<rsup|0>=p<rsup|0><rsub|max>-z \<xi\> <frac|q|2> ,
  </equation>

  while <math|p\<cdot\>k> is obtained from

  <\equation>
    (q-k<rsub|rec>)<rsup|2>=(p+k)<rsup|2>=m<rsup|2>+2p\<cdot\>k ,
  </equation>

  which yields

  <\eqnarray*>
    <tformat|<table|<row|<cell|q<rsup|2>+m<rsub|rec><rsup|2>-2q
    k<rsub|rec><rsup|0>>|<cell|=>|<cell|m<rsup|2>+2p\<cdot\>k
    ,>>|<row|<cell|p\<cdot\>k>|<cell|=>|<cell|<frac|q<rsup|2>+m<rsub|rec><rsup|2>-m<rsup|2>|2>-q
    k<rsub|rec> ,>>|<row|<cell|p\<cdot\>k>|<cell|=>|<cell|q(k<rsub|rec><rsup|max>-k<rsub|rec>)
    .<eq-number>>>>>
  </eqnarray*>

  Now

  <\equation>
    k<rsub|rec>=q-k<rsup|0>-p<rsup|0>=q-<left|(>\<xi\><frac|q|2>+p<rsub|max><rsup|0>-z\<xi\><frac|q|2><right|)>
    ,
  </equation>

  while <math|k<rsub|rec><rsup|max>> is obtained by setting <math|\<xi\>=0>
  in the above expression. Thus

  <\equation>
    p\<cdot\>k=q(k<rsub|rec><rsup|max>-k<rsub|rec>)=\<xi\>
    <frac|q<rsup|2>|2>(1-z) .
  </equation>

  Summarizing

  <\equation>
    K<rsub|\<perp\>><rsup|2>=<frac|\<xi\><rsup|2>
    q<rsup|3><left|(>1-z<right|)> |2p<rsup|0><rsub|max>-z \<xi\> q> \ .
  </equation>

  This expression seems to be still a bit complex. One could further simplify
  it by neglecting the term proportional to <math|\<xi\>> in the denominator.
  On the other hand, a theta function

  <\equation>
    \<theta\>(K<rsub|\<perp\>><rsup|2>-t)=\<theta\><left|(>\<xi\><rsup|2>q<rsup|3>-2
    t p<rsub|max><rsup|0>-(\<xi\>q<rsup|2>-t) z \<xi\>q<right|)>
  </equation>

  is easily solved in <math|z>, so that keeping the full expression is in
  fact a viable option.

  It is useful to have \ a view of the constant <math|K<rsub|T><rsup|2>>
  curves in the Dalitz plane.<\float|float|tbh>
    <\big-figure|<postscript|constantkperp.eps|*1/2|||||>>
      \;

      Lines of constant <math|K<rsub|\<perp\>><rsup|2>>.
    </big-figure>

    We see that these lines are all decreasing as a function of
    <math|k<rsup|0>>, and they cross the boundary of the Dalitz region just
    one. It will be convenient, in our case, to consider the extended region

    <\equation>
      z<rsub|2>\<less\>z\<less\>z<rsub|1> ,
      \ \ \ \ \ \ \ \ \ 0\<less\>\<xi\>\<less\>\<xi\><rsub|max>
      ,<label|eq:extendedregion>
    </equation>

    where

    <\equation>
      \<xi\><rsub|max>=1-<frac|(m+m<rsub|rec>)<rsup|2>|q<rsup|2>> .
    </equation>

    Phase space points generated outside the real Dalitz region will be
    easily rejected. Now, we can prove that, in the extended region
    (<reference|eq:extendedregion>), the lines of constant
    <math|K<rsub|\<perp\>><rsup|2>> are monotonically decreasing in <math|z>.
    In fact, solving for <math|\<xi\>> at fixed <math|z>, we find one
    positive solution

    <\equation>
      \<xi\>=<frac|<sqrt|K<rsub|\<perp\>><rsup|2><left|(>K<rsub|\<perp\>><rsup|2>z<rsup|2>+8p<rsup|0><rsub|max>q
      (1-z) <right|)>>-K<rsub|\<perp\>><rsup|2>z|2q<rsup|2>(1-z)>
      .<label|eq:xivsk2>
    </equation>

    It can be easily checked that the above function has positive derivative
    with respect to <math|z> for <math|K<rsub|\<perp\>><rsup|2>\<less\>2
    p<rsup|0><rsub|max> q>, and negative derivative for
    <math|K<rsub|\<perp\>><rsup|2>\<gtr\>2 p<rsup|0><rsub|max> q> . For the
    special value <math|K<rsub|\<perp\>><rsup|2>=2p<rsup|0><rsub|max> q>, the
    above equation yields <math|\<xi\>=><with|mode|math|2p<rsup|0><rsub|max>
    /q>, independent on <math|z>. On the other hand, this value is above
    <math|\<xi\><rsub|max>>,

    <\equation*>
      2 p<rsup|0><rsub|max> /q=<frac|q<rsup|2>+m<rsup|2>-M<rsub|rec><rsup|2>|q<rsup|2>>\<gtr\>1-<frac|(m+M<rsub|rec>)<rsup|2>|q<rsup|2>>
    </equation*>

    \ and therefore is never reached. Furthermore, eq.
    (<reference|eq:xivsk2>) is an increasing function of
    <math|K<rsub|\<perp\>><rsup|2>>, so, for
    <math|K<rsub|\<perp\>><rsup|2>\<gtr\>2 p<rsup|0><rsub|max> q>
    <math|\<xi\>\<gtr\>\<xi\><rsub|max>>. Thus, the minimum value of
    <math|\<xi\>> at given <math|K<rsub|\<perp\>><rsup|2>> in region
    (<reference|eq:extendedregion>) is

    <\equation>
      \<xi\><rsub|min>(K<rsub|\<perp\>><rsup|2>)=<frac|<sqrt|K<rsub|\<perp\>><rsup|2><left|(>K<rsub|\<perp\>><rsup|2>z<rsub|2><rsup|2>+8p<rsup|0><rsub|max>q
      (1-z<rsub|2>) <right|)>>-K<rsub|\<perp\>><rsup|2>z<rsub|2>|2q<rsup|2>(1-z<rsub|2>)>
      .
    </equation>

    If <math|\<xi\><rsub|min>(K<rsub|\<perp\>><rsup|2>)\<gtr\>\<xi\><rsub|max>>,
    this <math|K<rsub|\<perp\>>> value is forbidden. The maximum value of
    <math|K<rsub|\<perp\>>> is easily obtained, since it corresponds to
    <math|\<xi\>=\<xi\><rsub|max>> and <math|z=z<rsub|2>>

    <\equation*>
      t<rsub|max>\<equiv\>max(K<rsub|\<perp\>><rsup|2>)=<frac|\<xi\><rsub|max><rsup|2>
      q<rsup|3><left|(>1-z<rsub|2><right|)> |2p<rsup|0><rsub|max>-z<rsub|2>
      \<xi\><rsub|max> q> \ .
    </equation*>

    The integral of a phase space function with a
    <math|K<rsub|\<perp\>><rsup|2>> cut, assuming
    <math|t\<less\>t<rsub|max>>, is thus given by

    <\eqnarray*>
      <tformat|<table|<row|<cell|<big|int><rsub|0><rsup|\<xi\><rsub|max>>d\<xi\>
      <big|int><rsub|z<rsub|2>><rsup|z<rsub|1>>d z
      \ \ \ \ \<theta\><left|(>K<rsub|\<perp\>><rsup|2>-t<right|)>
      f(\<xi\>,z)>|<cell|=>|<cell|<big|int><rsub|\<xi\><rsub|min>(t)><rsup|\<xi\><rsub|max>>d\<xi\><big|int><rsup|min(z<rsub|1>,z<rsub|max>(t,\<xi\>))><rsub|z<rsub|2>>d
      z \ \ \ \ f(\<xi\>,z) ,<eq-number>>>>>
    </eqnarray*>

    where

    <\equation>
      z<rsub|max>(t,\<xi\>)=<frac|\<xi\><rsup|2>q<rsup|3>-2 t
      p<rsup|0><rsub|max>|\<xi\> q(\<xi\>q<rsup|2>- t)> .
    </equation>

    It is now convenient to further break this integral. Since
    <math|z<rsub|max>> is an increasing function of <math|\<xi\>>, there is
    going to be a value <math|\<xi\><rsub|1>> such that for
    <math|\<xi\>\<less\>\<xi\><rsub|1>> we have
    <math|z<rsub|max>(t,\<xi\>)\<less\>z<rsub|1>>. It is obviously given by

    <\equation>
      \<xi\><rsub|1>(t)=<frac|<sqrt|t<left|(>t
      z<rsub|1><rsup|2>+8p<rsup|0><rsub|max>q (1-z<rsub|1>) <right|)>>-t
      z<rsub|1>|2q<rsup|2>(1-z<rsub|1>)> .
    </equation>

    So, we further break the phase space as

    <\equation>
      <big|int><rsub|\<xi\><rsub|min>(t)><rsup|min(\<xi\><rsub|1>(t),\<xi\><rsub|max>)>d\<xi\><big|int><rsup|z<rsub|max>(t,\<xi\>))><rsub|z<rsub|2>>d
      z \ \ \ \ f(\<xi\>,z) +\<theta\>(\<xi\><rsub|max>-\<xi\><rsub|1>(t))<big|int><rsub|\<xi\><rsub|1>(t)><rsup|\<xi\><rsub|max>>d\<xi\><big|int><rsup|z<rsub|1>><rsub|z<rsub|2>>d
      z \ \ \ \ f(\<xi\>,z) .
    </equation>

    \;
  </float>

  <section|Approximation to <math|R/B>.>

  The iconal approximation to the amplitude yields

  <\equation>
    <with|math-font|cal|A><rsub|R>=<with|math-font|cal|A><rsub|B><left|(><frac|p<rsup|\<mu\>>|p\<cdot\>k>-<frac|r<rsup|\<mu\>>|r\<cdot\>k>
    \<ldots\><right|)> ,
  </equation>

  and squaring

  <\equation>
    <with|math-font|cal|A><rsub|R><rsup|2>=<with|math-font|cal|A><rsub|B><rsup|2><left|(>-<frac|m<rsup|2>|p\<cdot\>k>+<frac|2p\<cdot\>r|p\<cdot\>k
    \ \ \ \ \ r\<cdot\>k> \<ldots\><right|)>
  </equation>

  with the minus sign arising from the numerator of the gluon propagator. The
  <math|m<rsup|2>> term cannot \ make the cross section negative, and thus is
  bounded by the other terms. Separating the collinear region in the usual
  way

  <\equation*>
    <frac|p\<cdot\>r|p\<cdot\>k \ r\<cdot\>k>=<frac|p\<cdot\>r
    \ (p+r)\<cdot\>k|p\<cdot\>k \ \ r\<cdot\>k(p+r)\<cdot\>k>=<frac|p\<cdot\>r|p\<cdot\>k
    \ (p+r)\<cdot\>k>\<approx\><frac|E<rsup|2>|p\<cdot\>k E k<rsup|0>>
    =<frac|E|k<rsup|0> p\<cdot\>k>
  </equation*>

  <\equation>
    =<frac|1|\<xi\><rsup|2> q<rsup|2><left|(>1-z<right|)>><with|mode|text|><frac||><with|mode|text|>
    ,<label|eq:robapprox>
  </equation>

  where we drop constant numerical factors, since they enter in a
  redefinition of the upper bound norm. The jacobian for real radiation is

  <\equation>
    <frac|1|(2\<pi\>)<rsup|3>> <frac|q<rsup|0>|2<wide|k|\<bar\>><rsub|n>>
    \<times\><frac|q<rsup|0>|2>\<times\>k<rsup|0> d\<xi\> d z d\<phi\> .
  </equation>

  We can thus assume as an upper bound to the <math|R/B> expression the form
  (<reference|eq:robapprox>), extended to the whole region
  <math|0\<less\>\<xi\>\<less\>1> and <math|z<rsub|2>\<less\>z\<less\>z<rsub|1>>.
  A useful further restriction is to impose

  <\equation>
    K<rsub|\<perp\>><rsup|2>\<leqslant\>q<rsup|2>.
  </equation>

  This follows easily from

  <\equation>
    K<rsub|T><rsup|2>=<frac|k<rsup|0>|p<rsup|0>> 2 k\<cdot\>p
    \<leqslant\><frac|k<rsup|0>|p<rsup|0>> 4 k<rsup|0>
    p<rsup|0>=4(k<rsup|0>)<rsup|2>=q<rsup|2>\<xi\> <rsup|2>.
  </equation>

  Summarizing, our upper bounding function must have the form

  <\equation>
    <frac|q<rsup|0>|<wide|k|\<bar\>><rsub|n>> \ <frac|1|\<xi\>
    <left|(>1-z<right|)>><with|mode|text|><frac||><with|mode|text|> \ d
    \<xi\> d z d \<phi\>
  </equation>

  <section|Integration>

  We now evaluate the integral

  <\eqnarray*>
    <tformat|<table|<row|<cell|I(t)>|<cell|=>|<cell|<big|int><rsub|\<xi\><rsub|min>(t)><rsup|min(\<xi\><rsub|1>(t),\<xi\><rsub|max>)>d\<xi\><big|int><rsup|z<rsub|max>(t,\<xi\>)><rsub|z<rsub|2>>
    \ \ \ <frac|d z|\<xi\>(1-z)> +\<theta\>(\<xi\><rsub|max>-\<xi\><rsub|1>(t))<big|int><rsub|\<xi\><rsub|1>(t)><rsup|\<xi\><rsub|max>>d\<xi\><big|int><rsup|z<rsub|1>><rsub|z<rsub|2>>d
    z<frac|1|\<xi\>(1-z)> <eq-number>>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsub|\<xi\><rsub|min>(t)><rsup|min(\<xi\><rsub|1>(t),\<xi\><rsub|max>)><frac|d\<xi\>|\<xi\>>
    log<frac|1-z<rsub|2>|1-z<rsub|max>(t,\<xi\>)>
    +\<theta\>(\<xi\><rsub|max>-\<xi\><rsub|1>(t))
    \ log<frac|\<xi\><rsub|max>|\<xi\><rsub|1>(t)>
    log<frac|1-z<rsub|2>|1-z<rsub|1>> .<eq-number>>>>>
  </eqnarray*>

  In the first integral, the argument of the logarithm is

  <\equation>
    1-z<rsub|max>(t,\<xi\>)=1-<frac|\<xi\><rsup|2>q<rsup|3>-2 t
    p<rsup|0><rsub|max>|\<xi\> q(\<xi\>q<rsup|2>- t)>
    =<frac|\<xi\>q<rsup|3>-\<xi\>q t-\<xi\><rsup|2>q<rsup|3>+2t
    p<rsup|0><rsub|max>|\<xi\>q(\<xi\>q<rsup|2>-t)> = <frac|2t
    p<rsup|0><rsub|max>-\<xi\>q t|\<xi\>q(\<xi\>q<rsup|2>-t)> .
  </equation>

  We can break up the logarithm

  <\equation*>
    log<frac|(1-z<rsub|2>)\<xi\>q(\<xi\>q<rsup|2>-t)|t(2p<rsup|0><rsub|max>-\<xi\>
    q)>=log<left|[>(1-z<rsub|2>)<frac|q|t><right|]>+log\<xi\>+log(\<xi\>q<rsup|2>-t)-log(2p<rsup|0><rsub|max>-\<xi\>
    q)
  </equation*>

  Thus all integrals have the form

  <\equation>
    <big|int><frac|d\<xi\>|\<xi\>> log(a+b\<xi\>),
  </equation>

  with the argument of the logarithm being positive in the integration range.
  It can thus be expressed in terms of a dilogarithm

  <\equation>
    <big|int><frac|d\<xi\>|\<xi\>> log(a+b\<xi\>)=G(a,b,\<xi\>)+const. ,
    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 
  </equation>

  with

  <\equation>
    \ G(a,b,\<xi\>)\<equiv\>log(a+b\<xi\>)*log<left|(>1-<frac|a+b\<xi\>|a><right|)>+Li<rsub|2><left|(><frac|a+b\<xi\>|a><right|)>
    \ \ \ \ \ \ for a\<less\>0,
  </equation>

  <\equation>
    G(a,b,\<xi\>)\<equiv\>log<mid|\|><frac|b\<xi\>|a><mid|\|> log a -
    Li<rsub|2><left|[>-<frac|b\<xi\>|a><right|]>+<frac|\<pi\><rsup|6>|6>
    \ \ \ \ \ \ \ \ \ for a\<gtr\>0 .
  </equation>

  We thus get

  <\eqnarray*>
    <tformat|<table|<row|<cell|I(t)>|<cell|=>|<cell|<left|[>log \<xi\> log
    <left|[>(1-z<rsub|2>)<frac|q|t><right|]>+ <frac|1|2>log<rsup|2>
    \<xi\>+G(-t,q<rsup|2>,\<xi\>)-G(2p<rsup|0><rsub|max>,-q,\<xi\>)<right|]><rsup|min(\<xi\><rsub|1>(t),\<xi\><rsub|max>)><rsub|\<xi\><rsub|min>>>>|<row|<cell|>|<cell|+>|<cell|\<theta\>(\<xi\><rsub|max>-\<xi\><rsub|1>(t))
    \ log<frac|\<xi\><rsub|max>|\<xi\><rsub|1>(t)>
    log<frac|1-z<rsub|2>|1-z<rsub|1>> .<eq-number>>>>>
  </eqnarray*>

  <section|Generation of <math|z> and <math|\<xi\>> at fixed <math|t>.>

  We now see how <math|z> and <math|\<xi\>> can be generated once <math|t>
  has been found. They are distributed according to

  <\equation>
    <frac|d\<xi\> d z|\<xi\>(1-z)> \<delta\><left|(><frac|\<xi\><rsup|2>q<rsup|3>(1-z)|2p<rsup|0><rsub|max>-z\<xi\>q>-t<right|)>
    .
  </equation>

  Performing first the z integration, we get

  <\equation>
    <frac|d\<xi\>|\<xi\>><frac|1|(1-z)<mid|\|><frac|d|d z>
    <frac|\<xi\><rsup|2>q<rsup|3>(1-z)|2p<rsup|0><rsub|max>-z\<xi\>q><mid|\|>
    > ,
  </equation>

  evaluated at the <math|z> value that satisfies the <math|\<delta\>>
  function. We easily get

  <\equation>
    <frac|d\<xi\>|\<xi\>><big|int><frac| d z|(1-z)>
    \<delta\><left|(><frac|\<xi\><rsup|2>q<rsup|3>(1-z)|2p<rsup|0><rsub|max>-z\<xi\>q>-t<right|)>=
    \ d\<xi\> <frac|q<rsup|2>|t(\<xi\>q<rsup|2>-t)> .
  </equation>

  The <math|\<xi\>> value must thus be generated uniformly in
  <math|log(\<xi\>q<rsup|2>-t)>. The extremes for <math|\<xi\>> are given by
  <math|\<xi\><rsub|min>(t)> and <math|\<xi\><rsub|m>(t)=min(\<xi\><rsub|max>,\<xi\><rsub|1>(t))>

  <\equation>
    \<xi\>=<left|{>exp<left|[>log(\<xi\><rsub|min>(t) q<rsup|2>-t)+r
    log<frac|\<xi\><rsub|m>(t) q<rsup|2>-t|\<xi\><rsub|min>(t)q<rsup|2>-t><right|]>+t<right|}>/q<rsup|2>,
  </equation>

  where <math|0\<less\>r\<less\>1> is a uniform random number.

  <section|Summary>
</body>

<\initial>
  <\collection>
    <associate|sfactor|5>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|2>>
    <associate|auto-2|<tuple|1|3>>
    <associate|auto-3|<tuple|2|3>>
    <associate|auto-4|<tuple|2|5>>
    <associate|auto-5|<tuple|3|6>>
    <associate|auto-6|<tuple|4|7>>
    <associate|auto-7|<tuple|5|7>>
    <associate|auto-8|<tuple|6|?>>
    <associate|eq:extendedregion|<tuple|33|5>>
    <associate|eq:kfromy|<tuple|10|2>>
    <associate|eq:knp1max|<tuple|15|2>>
    <associate|eq:robapprox|<tuple|43|6>>
    <associate|eq:triangleeq|<tuple|11|2>>
    <associate|eq:xivsk2|<tuple|35|5>>
    <associate|fig:dalitz|<tuple|1|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<label|fig:dalitz>Dalitz region for
      <with|mode|<quote|math>|M<rsub|rec>= 0.3>,
      <with|mode|<quote|math>|m=0.2> .|<pageref|auto-2>>

      <\tuple|normal>
        \;

        Lines of constant <with|mode|<quote|math>|K<rsub|\<perp\>><rsup|2>>.
      </tuple|<pageref|auto-4>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Mapping
      to <with|mode|<quote|math>|><with|mode|<quote|math>|k<rsub|n+1><rsup|0>>
      to <with|mode|<quote|math>|y>> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Transverse
      momentum definition> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Approximation
      to <with|mode|<quote|math>|R/B>.> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Integration>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Generation
      of <with|mode|<quote|math>|z> and <with|mode|<quote|math>|\<xi\>> at
      fixed <with|mode|<quote|math>|t>.> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Summary>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>