      INTEGER process
      common /proc/ process

      DOUBLE PRECISION symmetryfactor
      common /factor/ symmetryfactor

      INTEGER narrowwidth
      common /narrow/ narrowwidth

      LOGICAL LOplusjet
      common /proc2/ LOplusjet

      INTEGER finallep, finallep_default, N_lepton
      common /leptonnumbers/ finallep(6), finallep_default(6), N_lepton

      INTEGER decaymode
      common /decaymode/ decaymode     
      
      INTEGER Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar
      INTEGER Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu
      INTEGER HjjjLO, HjjjLO_AA, HjjjLO_mu, HjjjLO_tau, HjjjLO_bbar
      INTEGER HjjjLO_WW, HjjjLO_ZZ_ll, HjjjLO_ZZ_lnu
      PARAMETER( Hjj           = 100  )
      PARAMETER( Hjj_AA        = 101  )
      PARAMETER( Hjj_mu        = 102  )
      PARAMETER( Hjj_tau       = 103  )
      PARAMETER( Hjj_bbar      = 104  )
      PARAMETER( Hjj_WW        = 105  )
      PARAMETER( Hjj_ZZ_ll     = 106  )
      PARAMETER( Hjj_ZZ_lnu    = 107  )
      PARAMETER( HjjjLO        = 9100 )
      PARAMETER( HjjjLO_AA     = 9101 )
      PARAMETER( HjjjLO_mu     = 9102 )
      PARAMETER( HjjjLO_tau    = 9103 )
      PARAMETER( HjjjLO_bbar   = 9104 )
      PARAMETER( HjjjLO_WW     = 9105 )
      PARAMETER( HjjjLO_ZZ_ll  = 9106 )
      PARAMETER( HjjjLO_ZZ_lnu = 9107 )

      INTEGER HJJJ,HJJJJLO
      INTEGER Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar
      INTEGER Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu
      INTEGER HjjjjLO_AA, HjjjjLO_mu, HjjjjLO_tau, HjjjjLO_bbar
      INTEGER HjjjjLO_WW, HjjjjLO_ZZ_ll, HjjjjLO_ZZ_lnu
      PARAMETER( Hjjj           = 110 )
      PARAMETER( Hjjj_AA        = 111  )
      PARAMETER( Hjjj_mu        = 112  )
      PARAMETER( Hjjj_tau       = 113  )
      PARAMETER( Hjjj_bbar      = 114  )
      PARAMETER( Hjjj_WW        = 115  )
      PARAMETER( Hjjj_ZZ_ll     = 116  )
      PARAMETER( Hjjj_ZZ_lnu    = 117  )
      PARAMETER( HjjjjLO        = 9110 )
      PARAMETER( HjjjjLO_AA     = 9111  )
      PARAMETER( HjjjjLO_mu     = 9112  )
      PARAMETER( HjjjjLO_tau    = 9113  )
      PARAMETER( HjjjjLO_bbar   = 9114  )
      PARAMETER( HjjjjLO_WW     = 9115  )
      PARAMETER( HjjjjLO_ZZ_ll  = 9116  )
      PARAMETER( HjjjjLO_ZZ_lnu = 9117  )

      INTEGER Zjj_l, Zjj_nu, WPjj, WMjj
      INTEGER ZjjjLO_l, ZjjjLO_nu, WPjjjLO, WMjjjLO
      INTEGER Ajj
      INTEGER AjjjLO
      PARAMETER( Zjj_l         = 120  )
      PARAMETER( Zjj_nu        = 121  )
      PARAMETER( WPjj          = 130  )
      PARAMETER( WMjj          = 140  )      
      PARAMETER( Ajj           = 150  )
      PARAMETER( ZjjjLO_l      = 9120 )
      PARAMETER( ZjjjLO_nu     = 9121 )
      PARAMETER( WPjjjLO       = 9130 )
      PARAMETER( WMjjjLO       = 9140 )
      PARAMETER( AjjjLO        = 9150 )

      INTEGER WPWMjj, WPWPjj, WMWMjj, WPZjj, WMZjj
      INTEGER ZZjj_ll, ZZjj_lnu,  AAjj
      INTEGER WPWMjjjLO,WPWPjjjLO,WMWMjjjLO,WPZjjjLO,WMZjjjLO
      INTEGER ZZjjjLO_ll,ZZjjjLO_lnu, AAjjjLO
      PARAMETER( WPWMjj        = 200  )
      PARAMETER( ZZjj_ll       = 210  )
      PARAMETER( ZZjj_lnu      = 211  )
      PARAMETER( WPZjj         = 220  )
      PARAMETER( WMZjj         = 230  )
      PARAMETER( AAjj          = 191  )
      PARAMETER( WPWPjj        = 250  )
      PARAMETER( WMWMjj        = 260  )
      PARAMETER( WPWMjjjLO     = 9200 )
      PARAMETER( ZZjjjLO_ll    = 9210 )
      PARAMETER( ZZjjjLO_lnu   = 9211 )
      PARAMETER( WPZjjjLO      = 9220 )
      PARAMETER( WMZjjjLO      = 9230 )
      PARAMETER( AAjjjLO       = 9191 )
      PARAMETER( WPWPjjjLO     = 9250  )
      PARAMETER( WMWMjjjLO     = 9260  )

      INTEGER WW, WPZ, WMZ, ZZ, WPA, WMA, ZA, AA,
     &        WWjLO, WPZjLOb, WMZjLOb, ZZjLO, WPAjLOb, WMAjLOb, ZAjLO, AAjLO
      PARAMETER( WW            = 300  )
      PARAMETER( WPZ           = 310  )
      PARAMETER( WMZ           = 320  )
      PARAMETER( ZZ            = 330  )
      PARAMETER( WPA           = 340  )
      PARAMETER( WMA           = 350  )
      PARAMETER( ZA            = 360  )
      PARAMETER( AA            = 370  )
      PARAMETER( WWjLO         = 9300 )
      PARAMETER( WPZjLOb       = 9310 )
      PARAMETER( WMZjLOb       = 9320 )
      PARAMETER( ZZjLO         = 9330 )
      PARAMETER( WPAjLOb       = 9340 )
      PARAMETER( WMAjLOb       = 9350 )       
      PARAMETER( ZAjLO         = 9360 )
      PARAMETER( AAjLO         = 9370 )

      INTEGER GFWW, GFZZ, GFZA, GFAA
      PARAMETER( GFWW          = 4300 )
      PARAMETER( GFZZ          = 4330 )
      PARAMETER( GFZA          = 4360 )
      PARAMETER( GFAA          = 4370 )


      INTEGER WWZ, ZZWP, ZZWM, WWWP, WWWM, WWA, ZZA, WPZA, WMZA,
     &        WPAA, WMAA, ZAA,ZnAA, AAA
      INTEGER WWZjLO, ZZWPjLO, ZZWMjLO, WWWPjLO, WWWMjLO, WWAjLO, ZZAjLO,
     &        WPZAjLO, WMZAjLO, WPAAjLO, WMAAjLO, ZAAjLO,ZnAAjLO, AAAjLO
      INTEGER WPAAj, WMAAj
      INTEGER WPAAjjLO, WMAAjjLO

      INTEGER ZZZ
      INTEGER ZZZjLO
      PARAMETER( WWZ           = 400  )
      PARAMETER( ZZWP          = 410  )
      PARAMETER( ZZWM          = 420  )
      PARAMETER( WWWP          = 430  )
      PARAMETER( WWWM          = 440  )
      PARAMETER( ZZZ           = 450  )
      PARAMETER( WWA           = 460  )
      PARAMETER( ZZA           = 470  )

      PARAMETER( WPZA          = 480  )
      PARAMETER( WMZA          = 490  )
      PARAMETER( WPAA          = 500  )
      PARAMETER( WMAA          = 510  )
      PARAMETER( ZAA           = 520  )
      PARAMETER( ZnAA          = 521  )
      PARAMETER( AAA           = 530  )
      PARAMETER( WWZjLO        = 9400 )
      PARAMETER( ZZWPjLO       = 9410 )
      PARAMETER( ZZWMjLO       = 9420 )
      PARAMETER( WWWPjLO       = 9430 )
      PARAMETER( WWWMjLO       = 9440 )
      PARAMETER( ZZZjLO        = 9450 )
      PARAMETER( WWAjLO        = 9460 )
      PARAMETER( ZZAjLO        = 9470 )

      PARAMETER( WPZAjLO       = 9480 )
      PARAMETER( WMZAjLO       = 9490 )
      PARAMETER( WPAAjLO       = 9500 )
      PARAMETER( WMAAjLO       = 9510 )
      PARAMETER( ZAAjLO        = 9520 )
      PARAMETER( ZnAAjLO       = 9521 )
      PARAMETER( AAAjLO        = 9530 )
      PARAMETER( WPAAj         = 800  )
      PARAMETER( WMAAj         = 810  )
      PARAMETER( WPAAjjLO      = 9800 )
      PARAMETER( WMAAjjLO      = 9810 )

c Higgs plus photon plus jets
      INTEGER HAjj, HAjj_AA, HAjj_mu, HAjj_tau, HAjj_bbar
      INTEGER HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu
      INTEGER HAjjjLO, HAjjjLO_AA, HAjjjLO_mu, HAjjjLO_tau, HAjjjLO_bbar
      INTEGER HAjjjLO_WW, HAjjjLO_ZZ_ll, HAjjjLO_ZZ_lnu

      PARAMETER( HAjj          = 2100 )
      PARAMETER( HAjj_AA       = 2101 )
      PARAMETER( HAjj_mu       = 2102 )
      PARAMETER( HAjj_tau      = 2103 )
      PARAMETER( HAjj_bbar     = 2104 )
      PARAMETER( HAjj_WW       = 2105 )
      PARAMETER( HAjj_ZZ_ll    = 2106 )
      PARAMETER( HAjj_ZZ_lnu   = 2107 )
      PARAMETER( HAjjjLO       = 92100)
      PARAMETER( HAjjjLO_AA    = 92101)
      PARAMETER( HAjjjLO_mu    = 92102)
      PARAMETER( HAjjjLO_tau   = 92103)
      PARAMETER( HAjjjLO_bbar  = 92104)
      PARAMETER( HAjjjLO_WW    = 92105)
      PARAMETER( HAjjjLO_ZZ_ll = 92106)
      PARAMETER( HAjjjLO_ZZ_lnu= 92107)
      
c diboson plus jet
      INTEGER WMZJ, WPZJ, WMAJ, WPAJ
      INTEGER WMZJLO, WPZJLO, WMAJLO, WPAJLO
      PARAMETER ( WMAJ   = 610 )
      PARAMETER ( WPAJ   = 620 )
      PARAMETER ( WMZJ   = 630 )
      PARAMETER ( WPZJ   = 640 )
      PARAMETER ( WMAJLO = 9610 )
      PARAMETER ( WPAJLO = 9620 )
      PARAMETER ( WMZJLO = 9630 )
      PARAMETER ( WPZJLO = 9640 )

c gluon fusion
      INTEGER GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar
      INTEGER GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu
      PARAMETER( GFHjj         = 4100 )
      PARAMETER( GFHjj_AA      = 4101 )
      PARAMETER( GFHjj_mu      = 4102 )
      PARAMETER( GFHjj_tau     = 4103 )
      PARAMETER( GFHjj_bbar    = 4104 )
      PARAMETER( GFHjj_WW      = 4105 )
      PARAMETER( GFHjj_ZZ_ll   = 4106 )
      PARAMETER( GFHjj_ZZ_lnu  = 4107 )

c -----------------------------------------
c processes not included in public version:
c -----------------------------------------

c QCD induced LO processes

      INTEGER QCDZjj_l, QCDZjj_nu, QCDWPWMjj, QCDZjjj, QCDWPjj
      INTEGER QCDWPjjj, QCDWMjj, QCDWMjjj
      PARAMETER( QCDZjj_l      = 5201 )
      PARAMETER( QCDZjj_nu     = 5204 )
      PARAMETER( QCDZjjj       = 5203 )
      PARAMETER( QCDWPjj       = 5302 )
      PARAMETER( QCDWPjjj      = 5303 )
      PARAMETER( QCDWMjj       = 5402 )
      PARAMETER( QCDWMjjj      = 5403 )
      PARAMETER( QCDWPWMjj     = 5411 )

      INTEGER EWtWj, EWtj, EWtZj_Wp, EWtZj_Wm
      PARAMETER( EWtWj         = 9414 )
      PARAMETER( EWtj          = 9415 )
      PARAMETER( EWtZj_Wp      = 9416 )
      PARAMETER( EWtZj_Wm      = 9417 )

c Drell-Yan and vector boson pair production

      INTEGER Z_DY, Wm_DY, Wp_DY
      INTEGER WPZkk, WMZkk
      PARAMETER(Z_DY           = 900  )
      PARAMETER(Wm_DY           = 910  )
      PARAMETER(Wp_DY           = 920  )  
      PARAMETER( WPZkk     = 311 )
      PARAMETER( WMZkk      = 321 )  
