      subroutine sreal_proc(p,legs,wgt)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision p(0:3,nexternal),wgt
      integer legs(nexternal),lstr
      character*140 str
      double precision P1(0:3,nexternal)
      integer i,ic(nexternal),legs1(nexternal)
      logical mtc,even
      
      do i=1,nexternal
         ic(i)=i
      enddo
      mtc=.false.
 10   call nexper(nexternal- 5,ic( 5+1),mtc,even)
      do i= 5+1,nexternal
         ic(i)=ic(i)+ 5
      enddo
      CALL SWITCHMOM (P,P1,IC,NEXTERNAL)
      CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL)
      
      call convert_to_string(nexternal,legs1,str,lstr)
      
      if (str.eq."-1-125-1112-1-2") then
         call srealmtrx_001(p1,wgt)
         goto 20
      elseif (str.eq."-1-125-1112-1-4") then
         call srealmtrx_002(p1,wgt)
         goto 20
      elseif (str.eq."-1125-11121-2") then
         call srealmtrx_003(p1,wgt)
         goto 20
      elseif (str.eq."-1125-11121-4") then
         call srealmtrx_004(p1,wgt)
         goto 20
      elseif (str.eq."-1125-1112-23") then
         call srealmtrx_005(p1,wgt)
         goto 20
      elseif (str.eq."-1125-1112-25") then
         call srealmtrx_006(p1,wgt)
         goto 20
      elseif (str.eq."-1125-1112-43") then
         call srealmtrx_007(p1,wgt)
         goto 20
      elseif (str.eq."-1125-1112-45") then
         call srealmtrx_008(p1,wgt)
         goto 20
      elseif (str.eq."-1-225-1112-2-2") then
         call srealmtrx_009(p1,wgt)
         goto 20
      elseif (str.eq."-1-225-1112-2-4") then
         call srealmtrx_010(p1,wgt)
         goto 20
      elseif (str.eq."-1225-11121-1") then
         call srealmtrx_011(p1,wgt)
         goto 20
      elseif (str.eq."-1225-1112-13") then
         call srealmtrx_012(p1,wgt)
         goto 20
      elseif (str.eq."-1225-1112-15") then
         call srealmtrx_013(p1,wgt)
         goto 20
      elseif (str.eq."-1225-11122-2") then
         call srealmtrx_014(p1,wgt)
         goto 20
      elseif (str.eq."-1225-11122-4") then
         call srealmtrx_015(p1,wgt)
         goto 20
      elseif (str.eq."-1225-11124-4") then
         call srealmtrx_016(p1,wgt)
         goto 20
      elseif (str.eq."-1225-11123-3") then
         call srealmtrx_017(p1,wgt)
         goto 20
      elseif (str.eq."-1225-11125-5") then
         call srealmtrx_018(p1,wgt)
         goto 20
      elseif (str.eq."-1225-111200") then
         call srealmtrx_019(p1,wgt)
         goto 20
      elseif (str.eq."-1-425-1112-2-4") then
         call srealmtrx_020(p1,wgt)
         goto 20
      elseif (str.eq."-1-425-1112-4-4") then
         call srealmtrx_021(p1,wgt)
         goto 20
      elseif (str.eq."-1425-11121-1") then
         call srealmtrx_022(p1,wgt)
         goto 20
      elseif (str.eq."-1425-1112-13") then
         call srealmtrx_023(p1,wgt)
         goto 20
      elseif (str.eq."-1425-1112-15") then
         call srealmtrx_024(p1,wgt)
         goto 20
      elseif (str.eq."-1425-11122-2") then
         call srealmtrx_025(p1,wgt)
         goto 20
      elseif (str.eq."-1425-1112-24") then
         call srealmtrx_026(p1,wgt)
         goto 20
      elseif (str.eq."-1425-11124-4") then
         call srealmtrx_027(p1,wgt)
         goto 20
      elseif (str.eq."-1425-11123-3") then
         call srealmtrx_028(p1,wgt)
         goto 20
      elseif (str.eq."-1425-11125-5") then
         call srealmtrx_029(p1,wgt)
         goto 20
      elseif (str.eq."-1425-111200") then
         call srealmtrx_030(p1,wgt)
         goto 20
      elseif (str.eq."-1-325-1112-1-2") then
         call srealmtrx_031(p1,wgt)
         goto 20
      elseif (str.eq."-1-325-1112-1-4") then
         call srealmtrx_032(p1,wgt)
         goto 20
      elseif (str.eq."-1-325-1112-2-3") then
         call srealmtrx_033(p1,wgt)
         goto 20
      elseif (str.eq."-1-325-1112-4-3") then
         call srealmtrx_034(p1,wgt)
         goto 20
      elseif (str.eq."-1325-1112-23") then
         call srealmtrx_035(p1,wgt)
         goto 20
      elseif (str.eq."-1325-1112-43") then
         call srealmtrx_036(p1,wgt)
         goto 20
      elseif (str.eq."-1-525-1112-1-2") then
         call srealmtrx_037(p1,wgt)
         goto 20
      elseif (str.eq."-1-525-1112-1-4") then
         call srealmtrx_038(p1,wgt)
         goto 20
      elseif (str.eq."-1-525-1112-2-5") then
         call srealmtrx_039(p1,wgt)
         goto 20
      elseif (str.eq."-1-525-1112-4-5") then
         call srealmtrx_040(p1,wgt)
         goto 20
      elseif (str.eq."-1525-1112-25") then
         call srealmtrx_041(p1,wgt)
         goto 20
      elseif (str.eq."-1525-1112-45") then
         call srealmtrx_042(p1,wgt)
         goto 20
      elseif (str.eq."-1025-1112-20") then
         call srealmtrx_043(p1,wgt)
         goto 20
      elseif (str.eq."-1025-1112-40") then
         call srealmtrx_044(p1,wgt)
         goto 20
      elseif (str.eq."1-125-11121-2") then
         call srealmtrx_045(p1,wgt)
         goto 20
      elseif (str.eq."1-125-11121-4") then
         call srealmtrx_046(p1,wgt)
         goto 20
      elseif (str.eq."1-125-1112-23") then
         call srealmtrx_047(p1,wgt)
         goto 20
      elseif (str.eq."1-125-1112-25") then
         call srealmtrx_048(p1,wgt)
         goto 20
      elseif (str.eq."1-125-1112-43") then
         call srealmtrx_049(p1,wgt)
         goto 20
      elseif (str.eq."1-125-1112-45") then
         call srealmtrx_050(p1,wgt)
         goto 20
      elseif (str.eq."1225-111211") then
         call srealmtrx_051(p1,wgt)
         goto 20
      elseif (str.eq."1225-111213") then
         call srealmtrx_052(p1,wgt)
         goto 20
      elseif (str.eq."1225-111215") then
         call srealmtrx_053(p1,wgt)
         goto 20
      elseif (str.eq."1425-111211") then
         call srealmtrx_054(p1,wgt)
         goto 20
      elseif (str.eq."1425-111213") then
         call srealmtrx_055(p1,wgt)
         goto 20
      elseif (str.eq."1425-111215") then
         call srealmtrx_056(p1,wgt)
         goto 20
      elseif (str.eq."1-325-11121-2") then
         call srealmtrx_057(p1,wgt)
         goto 20
      elseif (str.eq."1-325-11121-4") then
         call srealmtrx_058(p1,wgt)
         goto 20
      elseif (str.eq."1-525-11121-2") then
         call srealmtrx_059(p1,wgt)
         goto 20
      elseif (str.eq."1-525-11121-4") then
         call srealmtrx_060(p1,wgt)
         goto 20
      elseif (str.eq."-2-125-1112-2-2") then
         call srealmtrx_061(p1,wgt)
         goto 20
      elseif (str.eq."-2-125-1112-2-4") then
         call srealmtrx_062(p1,wgt)
         goto 20
      elseif (str.eq."-2225-11121-2") then
         call srealmtrx_063(p1,wgt)
         goto 20
      elseif (str.eq."-2225-11121-4") then
         call srealmtrx_064(p1,wgt)
         goto 20
      elseif (str.eq."-2225-1112-23") then
         call srealmtrx_065(p1,wgt)
         goto 20
      elseif (str.eq."-2225-1112-25") then
         call srealmtrx_066(p1,wgt)
         goto 20
      elseif (str.eq."-2225-1112-43") then
         call srealmtrx_067(p1,wgt)
         goto 20
      elseif (str.eq."-2225-1112-45") then
         call srealmtrx_068(p1,wgt)
         goto 20
      elseif (str.eq."-2425-11121-2") then
         call srealmtrx_069(p1,wgt)
         goto 20
      elseif (str.eq."-2425-1112-23") then
         call srealmtrx_070(p1,wgt)
         goto 20
      elseif (str.eq."-2425-1112-25") then
         call srealmtrx_071(p1,wgt)
         goto 20
      elseif (str.eq."-2-325-1112-2-2") then
         call srealmtrx_072(p1,wgt)
         goto 20
      elseif (str.eq."-2-325-1112-2-4") then
         call srealmtrx_073(p1,wgt)
         goto 20
      elseif (str.eq."-2-525-1112-2-2") then
         call srealmtrx_074(p1,wgt)
         goto 20
      elseif (str.eq."-2-525-1112-2-4") then
         call srealmtrx_075(p1,wgt)
         goto 20
      elseif (str.eq."2-125-11121-1") then
         call srealmtrx_076(p1,wgt)
         goto 20
      elseif (str.eq."2-125-1112-13") then
         call srealmtrx_077(p1,wgt)
         goto 20
      elseif (str.eq."2-125-1112-15") then
         call srealmtrx_078(p1,wgt)
         goto 20
      elseif (str.eq."2-125-11122-2") then
         call srealmtrx_079(p1,wgt)
         goto 20
      elseif (str.eq."2-125-11122-4") then
         call srealmtrx_080(p1,wgt)
         goto 20
      elseif (str.eq."2-125-11124-4") then
         call srealmtrx_081(p1,wgt)
         goto 20
      elseif (str.eq."2-125-11123-3") then
         call srealmtrx_082(p1,wgt)
         goto 20
      elseif (str.eq."2-125-11125-5") then
         call srealmtrx_083(p1,wgt)
         goto 20
      elseif (str.eq."2-125-111200") then
         call srealmtrx_084(p1,wgt)
         goto 20
      elseif (str.eq."2125-111211") then
         call srealmtrx_085(p1,wgt)
         goto 20
      elseif (str.eq."2125-111213") then
         call srealmtrx_086(p1,wgt)
         goto 20
      elseif (str.eq."2125-111215") then
         call srealmtrx_087(p1,wgt)
         goto 20
      elseif (str.eq."2-225-11121-2") then
         call srealmtrx_088(p1,wgt)
         goto 20
      elseif (str.eq."2-225-11121-4") then
         call srealmtrx_089(p1,wgt)
         goto 20
      elseif (str.eq."2-225-1112-23") then
         call srealmtrx_090(p1,wgt)
         goto 20
      elseif (str.eq."2-225-1112-25") then
         call srealmtrx_091(p1,wgt)
         goto 20
      elseif (str.eq."2-225-1112-43") then
         call srealmtrx_092(p1,wgt)
         goto 20
      elseif (str.eq."2-225-1112-45") then
         call srealmtrx_093(p1,wgt)
         goto 20
      elseif (str.eq."2225-111212") then
         call srealmtrx_094(p1,wgt)
         goto 20
      elseif (str.eq."2225-111223") then
         call srealmtrx_095(p1,wgt)
         goto 20
      elseif (str.eq."2225-111225") then
         call srealmtrx_096(p1,wgt)
         goto 20
      elseif (str.eq."2-425-11121-4") then
         call srealmtrx_097(p1,wgt)
         goto 20
      elseif (str.eq."2-425-1112-43") then
         call srealmtrx_098(p1,wgt)
         goto 20
      elseif (str.eq."2-425-1112-45") then
         call srealmtrx_099(p1,wgt)
         goto 20
      elseif (str.eq."2425-111212") then
         call srealmtrx_100(p1,wgt)
         goto 20
      elseif (str.eq."2425-111214") then
         call srealmtrx_101(p1,wgt)
         goto 20
      elseif (str.eq."2425-111223") then
         call srealmtrx_102(p1,wgt)
         goto 20
      elseif (str.eq."2425-111225") then
         call srealmtrx_103(p1,wgt)
         goto 20
      elseif (str.eq."2425-111243") then
         call srealmtrx_104(p1,wgt)
         goto 20
      elseif (str.eq."2425-111245") then
         call srealmtrx_105(p1,wgt)
         goto 20
      elseif (str.eq."2-325-11121-1") then
         call srealmtrx_106(p1,wgt)
         goto 20
      elseif (str.eq."2-325-11121-3") then
         call srealmtrx_107(p1,wgt)
         goto 20
      elseif (str.eq."2-325-11122-2") then
         call srealmtrx_108(p1,wgt)
         goto 20
      elseif (str.eq."2-325-11122-4") then
         call srealmtrx_109(p1,wgt)
         goto 20
      elseif (str.eq."2-325-11124-4") then
         call srealmtrx_110(p1,wgt)
         goto 20
      elseif (str.eq."2-325-11123-3") then
         call srealmtrx_111(p1,wgt)
         goto 20
      elseif (str.eq."2-325-1112-35") then
         call srealmtrx_112(p1,wgt)
         goto 20
      elseif (str.eq."2-325-11125-5") then
         call srealmtrx_113(p1,wgt)
         goto 20
      elseif (str.eq."2-325-111200") then
         call srealmtrx_114(p1,wgt)
         goto 20
      elseif (str.eq."2325-111213") then
         call srealmtrx_115(p1,wgt)
         goto 20
      elseif (str.eq."2325-111233") then
         call srealmtrx_116(p1,wgt)
         goto 20
      elseif (str.eq."2325-111235") then
         call srealmtrx_117(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11121-1") then
         call srealmtrx_118(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11121-5") then
         call srealmtrx_119(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11122-2") then
         call srealmtrx_120(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11122-4") then
         call srealmtrx_121(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11124-4") then
         call srealmtrx_122(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11123-3") then
         call srealmtrx_123(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11123-5") then
         call srealmtrx_124(p1,wgt)
         goto 20
      elseif (str.eq."2-525-11125-5") then
         call srealmtrx_125(p1,wgt)
         goto 20
      elseif (str.eq."2-525-111200") then
         call srealmtrx_126(p1,wgt)
         goto 20
      elseif (str.eq."2525-111215") then
         call srealmtrx_127(p1,wgt)
         goto 20
      elseif (str.eq."2525-111235") then
         call srealmtrx_128(p1,wgt)
         goto 20
      elseif (str.eq."2525-111255") then
         call srealmtrx_129(p1,wgt)
         goto 20
      elseif (str.eq."2025-111210") then
         call srealmtrx_130(p1,wgt)
         goto 20
      elseif (str.eq."2025-111230") then
         call srealmtrx_131(p1,wgt)
         goto 20
      elseif (str.eq."2025-111250") then
         call srealmtrx_132(p1,wgt)
         goto 20
      elseif (str.eq."-4-125-1112-2-4") then
         call srealmtrx_133(p1,wgt)
         goto 20
      elseif (str.eq."-4-125-1112-4-4") then
         call srealmtrx_134(p1,wgt)
         goto 20
      elseif (str.eq."-4225-11121-4") then
         call srealmtrx_135(p1,wgt)
         goto 20
      elseif (str.eq."-4225-1112-43") then
         call srealmtrx_136(p1,wgt)
         goto 20
      elseif (str.eq."-4225-1112-45") then
         call srealmtrx_137(p1,wgt)
         goto 20
      elseif (str.eq."-4425-11121-2") then
         call srealmtrx_138(p1,wgt)
         goto 20
      elseif (str.eq."-4425-11121-4") then
         call srealmtrx_139(p1,wgt)
         goto 20
      elseif (str.eq."-4425-1112-23") then
         call srealmtrx_140(p1,wgt)
         goto 20
      elseif (str.eq."-4425-1112-25") then
         call srealmtrx_141(p1,wgt)
         goto 20
      elseif (str.eq."-4425-1112-43") then
         call srealmtrx_142(p1,wgt)
         goto 20
      elseif (str.eq."-4425-1112-45") then
         call srealmtrx_143(p1,wgt)
         goto 20
      elseif (str.eq."-4-325-1112-2-4") then
         call srealmtrx_144(p1,wgt)
         goto 20
      elseif (str.eq."-4-325-1112-4-4") then
         call srealmtrx_145(p1,wgt)
         goto 20
      elseif (str.eq."-4-525-1112-2-4") then
         call srealmtrx_146(p1,wgt)
         goto 20
      elseif (str.eq."-4-525-1112-4-4") then
         call srealmtrx_147(p1,wgt)
         goto 20
      elseif (str.eq."4-125-11121-1") then
         call srealmtrx_148(p1,wgt)
         goto 20
      elseif (str.eq."4-125-1112-13") then
         call srealmtrx_149(p1,wgt)
         goto 20
      elseif (str.eq."4-125-1112-15") then
         call srealmtrx_150(p1,wgt)
         goto 20
      elseif (str.eq."4-125-11122-2") then
         call srealmtrx_151(p1,wgt)
         goto 20
      elseif (str.eq."4-125-1112-24") then
         call srealmtrx_152(p1,wgt)
         goto 20
      elseif (str.eq."4-125-11124-4") then
         call srealmtrx_153(p1,wgt)
         goto 20
      elseif (str.eq."4-125-11123-3") then
         call srealmtrx_154(p1,wgt)
         goto 20
      elseif (str.eq."4-125-11125-5") then
         call srealmtrx_155(p1,wgt)
         goto 20
      elseif (str.eq."4-125-111200") then
         call srealmtrx_156(p1,wgt)
         goto 20
      elseif (str.eq."4125-111211") then
         call srealmtrx_157(p1,wgt)
         goto 20
      elseif (str.eq."4125-111213") then
         call srealmtrx_158(p1,wgt)
         goto 20
      elseif (str.eq."4125-111215") then
         call srealmtrx_159(p1,wgt)
         goto 20
      elseif (str.eq."4-225-11121-2") then
         call srealmtrx_160(p1,wgt)
         goto 20
      elseif (str.eq."4-225-1112-23") then
         call srealmtrx_161(p1,wgt)
         goto 20
      elseif (str.eq."4-225-1112-25") then
         call srealmtrx_162(p1,wgt)
         goto 20
      elseif (str.eq."4225-111212") then
         call srealmtrx_163(p1,wgt)
         goto 20
      elseif (str.eq."4225-111214") then
         call srealmtrx_164(p1,wgt)
         goto 20
      elseif (str.eq."4225-111223") then
         call srealmtrx_165(p1,wgt)
         goto 20
      elseif (str.eq."4225-111225") then
         call srealmtrx_166(p1,wgt)
         goto 20
      elseif (str.eq."4225-111243") then
         call srealmtrx_167(p1,wgt)
         goto 20
      elseif (str.eq."4225-111245") then
         call srealmtrx_168(p1,wgt)
         goto 20
      elseif (str.eq."4-425-11121-2") then
         call srealmtrx_169(p1,wgt)
         goto 20
      elseif (str.eq."4-425-11121-4") then
         call srealmtrx_170(p1,wgt)
         goto 20
      elseif (str.eq."4-425-1112-23") then
         call srealmtrx_171(p1,wgt)
         goto 20
      elseif (str.eq."4-425-1112-25") then
         call srealmtrx_172(p1,wgt)
         goto 20
      elseif (str.eq."4-425-1112-43") then
         call srealmtrx_173(p1,wgt)
         goto 20
      elseif (str.eq."4-425-1112-45") then
         call srealmtrx_174(p1,wgt)
         goto 20
      elseif (str.eq."4425-111214") then
         call srealmtrx_175(p1,wgt)
         goto 20
      elseif (str.eq."4425-111243") then
         call srealmtrx_176(p1,wgt)
         goto 20
      elseif (str.eq."4425-111245") then
         call srealmtrx_177(p1,wgt)
         goto 20
      elseif (str.eq."4-325-11121-1") then
         call srealmtrx_178(p1,wgt)
         goto 20
      elseif (str.eq."4-325-11121-3") then
         call srealmtrx_179(p1,wgt)
         goto 20
      elseif (str.eq."4-325-11122-2") then
         call srealmtrx_180(p1,wgt)
         goto 20
      elseif (str.eq."4-325-1112-24") then
         call srealmtrx_181(p1,wgt)
         goto 20
      elseif (str.eq."4-325-11124-4") then
         call srealmtrx_182(p1,wgt)
         goto 20
      elseif (str.eq."4-325-11123-3") then
         call srealmtrx_183(p1,wgt)
         goto 20
      elseif (str.eq."4-325-1112-35") then
         call srealmtrx_184(p1,wgt)
         goto 20
      elseif (str.eq."4-325-11125-5") then
         call srealmtrx_185(p1,wgt)
         goto 20
      elseif (str.eq."4-325-111200") then
         call srealmtrx_186(p1,wgt)
         goto 20
      elseif (str.eq."4325-111213") then
         call srealmtrx_187(p1,wgt)
         goto 20
      elseif (str.eq."4325-111233") then
         call srealmtrx_188(p1,wgt)
         goto 20
      elseif (str.eq."4325-111235") then
         call srealmtrx_189(p1,wgt)
         goto 20
      elseif (str.eq."4-525-11121-1") then
         call srealmtrx_190(p1,wgt)
         goto 20
      elseif (str.eq."4-525-11121-5") then
         call srealmtrx_191(p1,wgt)
         goto 20
      elseif (str.eq."4-525-11122-2") then
         call srealmtrx_192(p1,wgt)
         goto 20
      elseif (str.eq."4-525-1112-24") then
         call srealmtrx_193(p1,wgt)
         goto 20
      elseif (str.eq."4-525-11124-4") then
         call srealmtrx_194(p1,wgt)
         goto 20
      elseif (str.eq."4-525-11123-3") then
         call srealmtrx_195(p1,wgt)
         goto 20
      elseif (str.eq."4-525-11123-5") then
         call srealmtrx_196(p1,wgt)
         goto 20
      elseif (str.eq."4-525-11125-5") then
         call srealmtrx_197(p1,wgt)
         goto 20
      elseif (str.eq."4-525-111200") then
         call srealmtrx_198(p1,wgt)
         goto 20
      elseif (str.eq."4525-111215") then
         call srealmtrx_199(p1,wgt)
         goto 20
      elseif (str.eq."4525-111235") then
         call srealmtrx_200(p1,wgt)
         goto 20
      elseif (str.eq."4525-111255") then
         call srealmtrx_201(p1,wgt)
         goto 20
      elseif (str.eq."4025-111210") then
         call srealmtrx_202(p1,wgt)
         goto 20
      elseif (str.eq."4025-111230") then
         call srealmtrx_203(p1,wgt)
         goto 20
      elseif (str.eq."4025-111250") then
         call srealmtrx_204(p1,wgt)
         goto 20
      elseif (str.eq."-3-125-1112-1-2") then
         call srealmtrx_205(p1,wgt)
         goto 20
      elseif (str.eq."-3-125-1112-1-4") then
         call srealmtrx_206(p1,wgt)
         goto 20
      elseif (str.eq."-3-125-1112-2-3") then
         call srealmtrx_207(p1,wgt)
         goto 20
      elseif (str.eq."-3-125-1112-4-3") then
         call srealmtrx_208(p1,wgt)
         goto 20
      elseif (str.eq."-3125-11121-2") then
         call srealmtrx_209(p1,wgt)
         goto 20
      elseif (str.eq."-3125-11121-4") then
         call srealmtrx_210(p1,wgt)
         goto 20
      elseif (str.eq."-3-225-1112-2-2") then
         call srealmtrx_211(p1,wgt)
         goto 20
      elseif (str.eq."-3-225-1112-2-4") then
         call srealmtrx_212(p1,wgt)
         goto 20
      elseif (str.eq."-3225-11121-1") then
         call srealmtrx_213(p1,wgt)
         goto 20
      elseif (str.eq."-3225-11121-3") then
         call srealmtrx_214(p1,wgt)
         goto 20
      elseif (str.eq."-3225-11122-2") then
         call srealmtrx_215(p1,wgt)
         goto 20
      elseif (str.eq."-3225-11122-4") then
         call srealmtrx_216(p1,wgt)
         goto 20
      elseif (str.eq."-3225-11124-4") then
         call srealmtrx_217(p1,wgt)
         goto 20
      elseif (str.eq."-3225-11123-3") then
         call srealmtrx_218(p1,wgt)
         goto 20
      elseif (str.eq."-3225-1112-35") then
         call srealmtrx_219(p1,wgt)
         goto 20
      elseif (str.eq."-3225-11125-5") then
         call srealmtrx_220(p1,wgt)
         goto 20
      elseif (str.eq."-3225-111200") then
         call srealmtrx_221(p1,wgt)
         goto 20
      elseif (str.eq."-3-425-1112-2-4") then
         call srealmtrx_222(p1,wgt)
         goto 20
      elseif (str.eq."-3-425-1112-4-4") then
         call srealmtrx_223(p1,wgt)
         goto 20
      elseif (str.eq."-3425-11121-1") then
         call srealmtrx_224(p1,wgt)
         goto 20
      elseif (str.eq."-3425-11121-3") then
         call srealmtrx_225(p1,wgt)
         goto 20
      elseif (str.eq."-3425-11122-2") then
         call srealmtrx_226(p1,wgt)
         goto 20
      elseif (str.eq."-3425-1112-24") then
         call srealmtrx_227(p1,wgt)
         goto 20
      elseif (str.eq."-3425-11124-4") then
         call srealmtrx_228(p1,wgt)
         goto 20
      elseif (str.eq."-3425-11123-3") then
         call srealmtrx_229(p1,wgt)
         goto 20
      elseif (str.eq."-3425-1112-35") then
         call srealmtrx_230(p1,wgt)
         goto 20
      elseif (str.eq."-3425-11125-5") then
         call srealmtrx_231(p1,wgt)
         goto 20
      elseif (str.eq."-3425-111200") then
         call srealmtrx_232(p1,wgt)
         goto 20
      elseif (str.eq."-3-325-1112-2-3") then
         call srealmtrx_233(p1,wgt)
         goto 20
      elseif (str.eq."-3-325-1112-4-3") then
         call srealmtrx_234(p1,wgt)
         goto 20
      elseif (str.eq."-3325-11121-2") then
         call srealmtrx_235(p1,wgt)
         goto 20
      elseif (str.eq."-3325-11121-4") then
         call srealmtrx_236(p1,wgt)
         goto 20
      elseif (str.eq."-3325-1112-23") then
         call srealmtrx_237(p1,wgt)
         goto 20
      elseif (str.eq."-3325-1112-25") then
         call srealmtrx_238(p1,wgt)
         goto 20
      elseif (str.eq."-3325-1112-43") then
         call srealmtrx_239(p1,wgt)
         goto 20
      elseif (str.eq."-3325-1112-45") then
         call srealmtrx_240(p1,wgt)
         goto 20
      elseif (str.eq."-3-525-1112-2-3") then
         call srealmtrx_241(p1,wgt)
         goto 20
      elseif (str.eq."-3-525-1112-2-5") then
         call srealmtrx_242(p1,wgt)
         goto 20
      elseif (str.eq."-3-525-1112-4-3") then
         call srealmtrx_243(p1,wgt)
         goto 20
      elseif (str.eq."-3-525-1112-4-5") then
         call srealmtrx_244(p1,wgt)
         goto 20
      elseif (str.eq."-3525-1112-25") then
         call srealmtrx_245(p1,wgt)
         goto 20
      elseif (str.eq."-3525-1112-45") then
         call srealmtrx_246(p1,wgt)
         goto 20
      elseif (str.eq."-3025-1112-20") then
         call srealmtrx_247(p1,wgt)
         goto 20
      elseif (str.eq."-3025-1112-40") then
         call srealmtrx_248(p1,wgt)
         goto 20
      elseif (str.eq."3-125-1112-23") then
         call srealmtrx_249(p1,wgt)
         goto 20
      elseif (str.eq."3-125-1112-43") then
         call srealmtrx_250(p1,wgt)
         goto 20
      elseif (str.eq."3225-111213") then
         call srealmtrx_251(p1,wgt)
         goto 20
      elseif (str.eq."3225-111233") then
         call srealmtrx_252(p1,wgt)
         goto 20
      elseif (str.eq."3225-111235") then
         call srealmtrx_253(p1,wgt)
         goto 20
      elseif (str.eq."3425-111213") then
         call srealmtrx_254(p1,wgt)
         goto 20
      elseif (str.eq."3425-111233") then
         call srealmtrx_255(p1,wgt)
         goto 20
      elseif (str.eq."3425-111235") then
         call srealmtrx_256(p1,wgt)
         goto 20
      elseif (str.eq."3-325-11121-2") then
         call srealmtrx_257(p1,wgt)
         goto 20
      elseif (str.eq."3-325-11121-4") then
         call srealmtrx_258(p1,wgt)
         goto 20
      elseif (str.eq."3-325-1112-23") then
         call srealmtrx_259(p1,wgt)
         goto 20
      elseif (str.eq."3-325-1112-25") then
         call srealmtrx_260(p1,wgt)
         goto 20
      elseif (str.eq."3-325-1112-43") then
         call srealmtrx_261(p1,wgt)
         goto 20
      elseif (str.eq."3-325-1112-45") then
         call srealmtrx_262(p1,wgt)
         goto 20
      elseif (str.eq."3-525-1112-23") then
         call srealmtrx_263(p1,wgt)
         goto 20
      elseif (str.eq."3-525-1112-43") then
         call srealmtrx_264(p1,wgt)
         goto 20
      elseif (str.eq."-5-125-1112-1-2") then
         call srealmtrx_265(p1,wgt)
         goto 20
      elseif (str.eq."-5-125-1112-1-4") then
         call srealmtrx_266(p1,wgt)
         goto 20
      elseif (str.eq."-5-125-1112-2-5") then
         call srealmtrx_267(p1,wgt)
         goto 20
      elseif (str.eq."-5-125-1112-4-5") then
         call srealmtrx_268(p1,wgt)
         goto 20
      elseif (str.eq."-5125-11121-2") then
         call srealmtrx_269(p1,wgt)
         goto 20
      elseif (str.eq."-5125-11121-4") then
         call srealmtrx_270(p1,wgt)
         goto 20
      elseif (str.eq."-5-225-1112-2-2") then
         call srealmtrx_271(p1,wgt)
         goto 20
      elseif (str.eq."-5-225-1112-2-4") then
         call srealmtrx_272(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11121-1") then
         call srealmtrx_273(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11121-5") then
         call srealmtrx_274(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11122-2") then
         call srealmtrx_275(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11122-4") then
         call srealmtrx_276(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11124-4") then
         call srealmtrx_277(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11123-3") then
         call srealmtrx_278(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11123-5") then
         call srealmtrx_279(p1,wgt)
         goto 20
      elseif (str.eq."-5225-11125-5") then
         call srealmtrx_280(p1,wgt)
         goto 20
      elseif (str.eq."-5225-111200") then
         call srealmtrx_281(p1,wgt)
         goto 20
      elseif (str.eq."-5-425-1112-2-4") then
         call srealmtrx_282(p1,wgt)
         goto 20
      elseif (str.eq."-5-425-1112-4-4") then
         call srealmtrx_283(p1,wgt)
         goto 20
      elseif (str.eq."-5425-11121-1") then
         call srealmtrx_284(p1,wgt)
         goto 20
      elseif (str.eq."-5425-11121-5") then
         call srealmtrx_285(p1,wgt)
         goto 20
      elseif (str.eq."-5425-11122-2") then
         call srealmtrx_286(p1,wgt)
         goto 20
      elseif (str.eq."-5425-1112-24") then
         call srealmtrx_287(p1,wgt)
         goto 20
      elseif (str.eq."-5425-11124-4") then
         call srealmtrx_288(p1,wgt)
         goto 20
      elseif (str.eq."-5425-11123-3") then
         call srealmtrx_289(p1,wgt)
         goto 20
      elseif (str.eq."-5425-11123-5") then
         call srealmtrx_290(p1,wgt)
         goto 20
      elseif (str.eq."-5425-11125-5") then
         call srealmtrx_291(p1,wgt)
         goto 20
      elseif (str.eq."-5425-111200") then
         call srealmtrx_292(p1,wgt)
         goto 20
      elseif (str.eq."-5-325-1112-2-3") then
         call srealmtrx_293(p1,wgt)
         goto 20
      elseif (str.eq."-5-325-1112-2-5") then
         call srealmtrx_294(p1,wgt)
         goto 20
      elseif (str.eq."-5-325-1112-4-3") then
         call srealmtrx_295(p1,wgt)
         goto 20
      elseif (str.eq."-5-325-1112-4-5") then
         call srealmtrx_296(p1,wgt)
         goto 20
      elseif (str.eq."-5325-1112-23") then
         call srealmtrx_297(p1,wgt)
         goto 20
      elseif (str.eq."-5325-1112-43") then
         call srealmtrx_298(p1,wgt)
         goto 20
      elseif (str.eq."-5-525-1112-2-5") then
         call srealmtrx_299(p1,wgt)
         goto 20
      elseif (str.eq."-5-525-1112-4-5") then
         call srealmtrx_300(p1,wgt)
         goto 20
      elseif (str.eq."-5525-11121-2") then
         call srealmtrx_301(p1,wgt)
         goto 20
      elseif (str.eq."-5525-11121-4") then
         call srealmtrx_302(p1,wgt)
         goto 20
      elseif (str.eq."-5525-1112-23") then
         call srealmtrx_303(p1,wgt)
         goto 20
      elseif (str.eq."-5525-1112-25") then
         call srealmtrx_304(p1,wgt)
         goto 20
      elseif (str.eq."-5525-1112-43") then
         call srealmtrx_305(p1,wgt)
         goto 20
      elseif (str.eq."-5525-1112-45") then
         call srealmtrx_306(p1,wgt)
         goto 20
      elseif (str.eq."-5025-1112-20") then
         call srealmtrx_307(p1,wgt)
         goto 20
      elseif (str.eq."-5025-1112-40") then
         call srealmtrx_308(p1,wgt)
         goto 20
      elseif (str.eq."5-125-1112-25") then
         call srealmtrx_309(p1,wgt)
         goto 20
      elseif (str.eq."5-125-1112-45") then
         call srealmtrx_310(p1,wgt)
         goto 20
      elseif (str.eq."5225-111215") then
         call srealmtrx_311(p1,wgt)
         goto 20
      elseif (str.eq."5225-111235") then
         call srealmtrx_312(p1,wgt)
         goto 20
      elseif (str.eq."5225-111255") then
         call srealmtrx_313(p1,wgt)
         goto 20
      elseif (str.eq."5425-111215") then
         call srealmtrx_314(p1,wgt)
         goto 20
      elseif (str.eq."5425-111235") then
         call srealmtrx_315(p1,wgt)
         goto 20
      elseif (str.eq."5425-111255") then
         call srealmtrx_316(p1,wgt)
         goto 20
      elseif (str.eq."5-325-1112-25") then
         call srealmtrx_317(p1,wgt)
         goto 20
      elseif (str.eq."5-325-1112-45") then
         call srealmtrx_318(p1,wgt)
         goto 20
      elseif (str.eq."5-525-11121-2") then
         call srealmtrx_319(p1,wgt)
         goto 20
      elseif (str.eq."5-525-11121-4") then
         call srealmtrx_320(p1,wgt)
         goto 20
      elseif (str.eq."5-525-1112-23") then
         call srealmtrx_321(p1,wgt)
         goto 20
      elseif (str.eq."5-525-1112-25") then
         call srealmtrx_322(p1,wgt)
         goto 20
      elseif (str.eq."5-525-1112-43") then
         call srealmtrx_323(p1,wgt)
         goto 20
      elseif (str.eq."5-525-1112-45") then
         call srealmtrx_324(p1,wgt)
         goto 20
      elseif (str.eq."0-125-1112-20") then
         call srealmtrx_325(p1,wgt)
         goto 20
      elseif (str.eq."0-125-1112-40") then
         call srealmtrx_326(p1,wgt)
         goto 20
      elseif (str.eq."0225-111210") then
         call srealmtrx_327(p1,wgt)
         goto 20
      elseif (str.eq."0225-111230") then
         call srealmtrx_328(p1,wgt)
         goto 20
      elseif (str.eq."0225-111250") then
         call srealmtrx_329(p1,wgt)
         goto 20
      elseif (str.eq."0425-111210") then
         call srealmtrx_330(p1,wgt)
         goto 20
      elseif (str.eq."0425-111230") then
         call srealmtrx_331(p1,wgt)
         goto 20
      elseif (str.eq."0425-111250") then
         call srealmtrx_332(p1,wgt)
         goto 20
      elseif (str.eq."0-325-1112-20") then
         call srealmtrx_333(p1,wgt)
         goto 20
      elseif (str.eq."0-325-1112-40") then
         call srealmtrx_334(p1,wgt)
         goto 20
      elseif (str.eq."0-525-1112-20") then
         call srealmtrx_335(p1,wgt)
         goto 20
      elseif (str.eq."0-525-1112-40") then
         call srealmtrx_336(p1,wgt)
         goto 20
      elseif (str.eq."0025-11121-2") then
         call srealmtrx_337(p1,wgt)
         goto 20
      elseif (str.eq."0025-11121-4") then
         call srealmtrx_338(p1,wgt)
         goto 20
      elseif (str.eq."0025-1112-23") then
         call srealmtrx_339(p1,wgt)
         goto 20
      elseif (str.eq."0025-1112-25") then
         call srealmtrx_340(p1,wgt)
         goto 20
      elseif (str.eq."0025-1112-43") then
         call srealmtrx_341(p1,wgt)
         goto 20
      elseif (str.eq."0025-1112-45") then
         call srealmtrx_342(p1,wgt)
         goto 20
      endif
      
      do while(mtc)
         do i= 5+1,nexternal
            ic(i)=ic(i)- 5
         enddo
         goto 10
      enddo
      if (.not.mtc) then
         write (*,*) "Error #1, in sreal_proc.f"
         stop
      endif
      
 20   continue
      return
      end
      
      
      subroutine real_color(legs,color)
      implicit none
      include "nexternal.inc"
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2001(maxamps), jamp2001(0:maxamps)
      common/to_Ramps_001/amp2001,jamp2001
      Double Precision amp2002(maxamps), jamp2002(0:maxamps)
      common/to_Ramps_002/amp2002,jamp2002
      Double Precision amp2003(maxamps), jamp2003(0:maxamps)
      common/to_Ramps_003/amp2003,jamp2003
      Double Precision amp2004(maxamps), jamp2004(0:maxamps)
      common/to_Ramps_004/amp2004,jamp2004
      Double Precision amp2005(maxamps), jamp2005(0:maxamps)
      common/to_Ramps_005/amp2005,jamp2005
      Double Precision amp2006(maxamps), jamp2006(0:maxamps)
      common/to_Ramps_006/amp2006,jamp2006
      Double Precision amp2007(maxamps), jamp2007(0:maxamps)
      common/to_Ramps_007/amp2007,jamp2007
      Double Precision amp2008(maxamps), jamp2008(0:maxamps)
      common/to_Ramps_008/amp2008,jamp2008
      Double Precision amp2009(maxamps), jamp2009(0:maxamps)
      common/to_Ramps_009/amp2009,jamp2009
      Double Precision amp2010(maxamps), jamp2010(0:maxamps)
      common/to_Ramps_010/amp2010,jamp2010
      Double Precision amp2011(maxamps), jamp2011(0:maxamps)
      common/to_Ramps_011/amp2011,jamp2011
      Double Precision amp2012(maxamps), jamp2012(0:maxamps)
      common/to_Ramps_012/amp2012,jamp2012
      Double Precision amp2013(maxamps), jamp2013(0:maxamps)
      common/to_Ramps_013/amp2013,jamp2013
      Double Precision amp2014(maxamps), jamp2014(0:maxamps)
      common/to_Ramps_014/amp2014,jamp2014
      Double Precision amp2015(maxamps), jamp2015(0:maxamps)
      common/to_Ramps_015/amp2015,jamp2015
      Double Precision amp2016(maxamps), jamp2016(0:maxamps)
      common/to_Ramps_016/amp2016,jamp2016
      Double Precision amp2017(maxamps), jamp2017(0:maxamps)
      common/to_Ramps_017/amp2017,jamp2017
      Double Precision amp2018(maxamps), jamp2018(0:maxamps)
      common/to_Ramps_018/amp2018,jamp2018
      Double Precision amp2019(maxamps), jamp2019(0:maxamps)
      common/to_Ramps_019/amp2019,jamp2019
      Double Precision amp2020(maxamps), jamp2020(0:maxamps)
      common/to_Ramps_020/amp2020,jamp2020
      Double Precision amp2021(maxamps), jamp2021(0:maxamps)
      common/to_Ramps_021/amp2021,jamp2021
      Double Precision amp2022(maxamps), jamp2022(0:maxamps)
      common/to_Ramps_022/amp2022,jamp2022
      Double Precision amp2023(maxamps), jamp2023(0:maxamps)
      common/to_Ramps_023/amp2023,jamp2023
      Double Precision amp2024(maxamps), jamp2024(0:maxamps)
      common/to_Ramps_024/amp2024,jamp2024
      Double Precision amp2025(maxamps), jamp2025(0:maxamps)
      common/to_Ramps_025/amp2025,jamp2025
      Double Precision amp2026(maxamps), jamp2026(0:maxamps)
      common/to_Ramps_026/amp2026,jamp2026
      Double Precision amp2027(maxamps), jamp2027(0:maxamps)
      common/to_Ramps_027/amp2027,jamp2027
      Double Precision amp2028(maxamps), jamp2028(0:maxamps)
      common/to_Ramps_028/amp2028,jamp2028
      Double Precision amp2029(maxamps), jamp2029(0:maxamps)
      common/to_Ramps_029/amp2029,jamp2029
      Double Precision amp2030(maxamps), jamp2030(0:maxamps)
      common/to_Ramps_030/amp2030,jamp2030
      Double Precision amp2031(maxamps), jamp2031(0:maxamps)
      common/to_Ramps_031/amp2031,jamp2031
      Double Precision amp2032(maxamps), jamp2032(0:maxamps)
      common/to_Ramps_032/amp2032,jamp2032
      Double Precision amp2033(maxamps), jamp2033(0:maxamps)
      common/to_Ramps_033/amp2033,jamp2033
      Double Precision amp2034(maxamps), jamp2034(0:maxamps)
      common/to_Ramps_034/amp2034,jamp2034
      Double Precision amp2035(maxamps), jamp2035(0:maxamps)
      common/to_Ramps_035/amp2035,jamp2035
      Double Precision amp2036(maxamps), jamp2036(0:maxamps)
      common/to_Ramps_036/amp2036,jamp2036
      Double Precision amp2037(maxamps), jamp2037(0:maxamps)
      common/to_Ramps_037/amp2037,jamp2037
      Double Precision amp2038(maxamps), jamp2038(0:maxamps)
      common/to_Ramps_038/amp2038,jamp2038
      Double Precision amp2039(maxamps), jamp2039(0:maxamps)
      common/to_Ramps_039/amp2039,jamp2039
      Double Precision amp2040(maxamps), jamp2040(0:maxamps)
      common/to_Ramps_040/amp2040,jamp2040
      Double Precision amp2041(maxamps), jamp2041(0:maxamps)
      common/to_Ramps_041/amp2041,jamp2041
      Double Precision amp2042(maxamps), jamp2042(0:maxamps)
      common/to_Ramps_042/amp2042,jamp2042
      Double Precision amp2043(maxamps), jamp2043(0:maxamps)
      common/to_Ramps_043/amp2043,jamp2043
      Double Precision amp2044(maxamps), jamp2044(0:maxamps)
      common/to_Ramps_044/amp2044,jamp2044
      Double Precision amp2045(maxamps), jamp2045(0:maxamps)
      common/to_Ramps_045/amp2045,jamp2045
      Double Precision amp2046(maxamps), jamp2046(0:maxamps)
      common/to_Ramps_046/amp2046,jamp2046
      Double Precision amp2047(maxamps), jamp2047(0:maxamps)
      common/to_Ramps_047/amp2047,jamp2047
      Double Precision amp2048(maxamps), jamp2048(0:maxamps)
      common/to_Ramps_048/amp2048,jamp2048
      Double Precision amp2049(maxamps), jamp2049(0:maxamps)
      common/to_Ramps_049/amp2049,jamp2049
      Double Precision amp2050(maxamps), jamp2050(0:maxamps)
      common/to_Ramps_050/amp2050,jamp2050
      Double Precision amp2051(maxamps), jamp2051(0:maxamps)
      common/to_Ramps_051/amp2051,jamp2051
      Double Precision amp2052(maxamps), jamp2052(0:maxamps)
      common/to_Ramps_052/amp2052,jamp2052
      Double Precision amp2053(maxamps), jamp2053(0:maxamps)
      common/to_Ramps_053/amp2053,jamp2053
      Double Precision amp2054(maxamps), jamp2054(0:maxamps)
      common/to_Ramps_054/amp2054,jamp2054
      Double Precision amp2055(maxamps), jamp2055(0:maxamps)
      common/to_Ramps_055/amp2055,jamp2055
      Double Precision amp2056(maxamps), jamp2056(0:maxamps)
      common/to_Ramps_056/amp2056,jamp2056
      Double Precision amp2057(maxamps), jamp2057(0:maxamps)
      common/to_Ramps_057/amp2057,jamp2057
      Double Precision amp2058(maxamps), jamp2058(0:maxamps)
      common/to_Ramps_058/amp2058,jamp2058
      Double Precision amp2059(maxamps), jamp2059(0:maxamps)
      common/to_Ramps_059/amp2059,jamp2059
      Double Precision amp2060(maxamps), jamp2060(0:maxamps)
      common/to_Ramps_060/amp2060,jamp2060
      Double Precision amp2061(maxamps), jamp2061(0:maxamps)
      common/to_Ramps_061/amp2061,jamp2061
      Double Precision amp2062(maxamps), jamp2062(0:maxamps)
      common/to_Ramps_062/amp2062,jamp2062
      Double Precision amp2063(maxamps), jamp2063(0:maxamps)
      common/to_Ramps_063/amp2063,jamp2063
      Double Precision amp2064(maxamps), jamp2064(0:maxamps)
      common/to_Ramps_064/amp2064,jamp2064
      Double Precision amp2065(maxamps), jamp2065(0:maxamps)
      common/to_Ramps_065/amp2065,jamp2065
      Double Precision amp2066(maxamps), jamp2066(0:maxamps)
      common/to_Ramps_066/amp2066,jamp2066
      Double Precision amp2067(maxamps), jamp2067(0:maxamps)
      common/to_Ramps_067/amp2067,jamp2067
      Double Precision amp2068(maxamps), jamp2068(0:maxamps)
      common/to_Ramps_068/amp2068,jamp2068
      Double Precision amp2069(maxamps), jamp2069(0:maxamps)
      common/to_Ramps_069/amp2069,jamp2069
      Double Precision amp2070(maxamps), jamp2070(0:maxamps)
      common/to_Ramps_070/amp2070,jamp2070
      Double Precision amp2071(maxamps), jamp2071(0:maxamps)
      common/to_Ramps_071/amp2071,jamp2071
      Double Precision amp2072(maxamps), jamp2072(0:maxamps)
      common/to_Ramps_072/amp2072,jamp2072
      Double Precision amp2073(maxamps), jamp2073(0:maxamps)
      common/to_Ramps_073/amp2073,jamp2073
      Double Precision amp2074(maxamps), jamp2074(0:maxamps)
      common/to_Ramps_074/amp2074,jamp2074
      Double Precision amp2075(maxamps), jamp2075(0:maxamps)
      common/to_Ramps_075/amp2075,jamp2075
      Double Precision amp2076(maxamps), jamp2076(0:maxamps)
      common/to_Ramps_076/amp2076,jamp2076
      Double Precision amp2077(maxamps), jamp2077(0:maxamps)
      common/to_Ramps_077/amp2077,jamp2077
      Double Precision amp2078(maxamps), jamp2078(0:maxamps)
      common/to_Ramps_078/amp2078,jamp2078
      Double Precision amp2079(maxamps), jamp2079(0:maxamps)
      common/to_Ramps_079/amp2079,jamp2079
      Double Precision amp2080(maxamps), jamp2080(0:maxamps)
      common/to_Ramps_080/amp2080,jamp2080
      Double Precision amp2081(maxamps), jamp2081(0:maxamps)
      common/to_Ramps_081/amp2081,jamp2081
      Double Precision amp2082(maxamps), jamp2082(0:maxamps)
      common/to_Ramps_082/amp2082,jamp2082
      Double Precision amp2083(maxamps), jamp2083(0:maxamps)
      common/to_Ramps_083/amp2083,jamp2083
      Double Precision amp2084(maxamps), jamp2084(0:maxamps)
      common/to_Ramps_084/amp2084,jamp2084
      Double Precision amp2085(maxamps), jamp2085(0:maxamps)
      common/to_Ramps_085/amp2085,jamp2085
      Double Precision amp2086(maxamps), jamp2086(0:maxamps)
      common/to_Ramps_086/amp2086,jamp2086
      Double Precision amp2087(maxamps), jamp2087(0:maxamps)
      common/to_Ramps_087/amp2087,jamp2087
      Double Precision amp2088(maxamps), jamp2088(0:maxamps)
      common/to_Ramps_088/amp2088,jamp2088
      Double Precision amp2089(maxamps), jamp2089(0:maxamps)
      common/to_Ramps_089/amp2089,jamp2089
      Double Precision amp2090(maxamps), jamp2090(0:maxamps)
      common/to_Ramps_090/amp2090,jamp2090
      Double Precision amp2091(maxamps), jamp2091(0:maxamps)
      common/to_Ramps_091/amp2091,jamp2091
      Double Precision amp2092(maxamps), jamp2092(0:maxamps)
      common/to_Ramps_092/amp2092,jamp2092
      Double Precision amp2093(maxamps), jamp2093(0:maxamps)
      common/to_Ramps_093/amp2093,jamp2093
      Double Precision amp2094(maxamps), jamp2094(0:maxamps)
      common/to_Ramps_094/amp2094,jamp2094
      Double Precision amp2095(maxamps), jamp2095(0:maxamps)
      common/to_Ramps_095/amp2095,jamp2095
      Double Precision amp2096(maxamps), jamp2096(0:maxamps)
      common/to_Ramps_096/amp2096,jamp2096
      Double Precision amp2097(maxamps), jamp2097(0:maxamps)
      common/to_Ramps_097/amp2097,jamp2097
      Double Precision amp2098(maxamps), jamp2098(0:maxamps)
      common/to_Ramps_098/amp2098,jamp2098
      Double Precision amp2099(maxamps), jamp2099(0:maxamps)
      common/to_Ramps_099/amp2099,jamp2099
      Double Precision amp2100(maxamps), jamp2100(0:maxamps)
      common/to_Ramps_100/amp2100,jamp2100
      Double Precision amp2101(maxamps), jamp2101(0:maxamps)
      common/to_Ramps_101/amp2101,jamp2101
      Double Precision amp2102(maxamps), jamp2102(0:maxamps)
      common/to_Ramps_102/amp2102,jamp2102
      Double Precision amp2103(maxamps), jamp2103(0:maxamps)
      common/to_Ramps_103/amp2103,jamp2103
      Double Precision amp2104(maxamps), jamp2104(0:maxamps)
      common/to_Ramps_104/amp2104,jamp2104
      Double Precision amp2105(maxamps), jamp2105(0:maxamps)
      common/to_Ramps_105/amp2105,jamp2105
      Double Precision amp2106(maxamps), jamp2106(0:maxamps)
      common/to_Ramps_106/amp2106,jamp2106
      Double Precision amp2107(maxamps), jamp2107(0:maxamps)
      common/to_Ramps_107/amp2107,jamp2107
      Double Precision amp2108(maxamps), jamp2108(0:maxamps)
      common/to_Ramps_108/amp2108,jamp2108
      Double Precision amp2109(maxamps), jamp2109(0:maxamps)
      common/to_Ramps_109/amp2109,jamp2109
      Double Precision amp2110(maxamps), jamp2110(0:maxamps)
      common/to_Ramps_110/amp2110,jamp2110
      Double Precision amp2111(maxamps), jamp2111(0:maxamps)
      common/to_Ramps_111/amp2111,jamp2111
      Double Precision amp2112(maxamps), jamp2112(0:maxamps)
      common/to_Ramps_112/amp2112,jamp2112
      Double Precision amp2113(maxamps), jamp2113(0:maxamps)
      common/to_Ramps_113/amp2113,jamp2113
      Double Precision amp2114(maxamps), jamp2114(0:maxamps)
      common/to_Ramps_114/amp2114,jamp2114
      Double Precision amp2115(maxamps), jamp2115(0:maxamps)
      common/to_Ramps_115/amp2115,jamp2115
      Double Precision amp2116(maxamps), jamp2116(0:maxamps)
      common/to_Ramps_116/amp2116,jamp2116
      Double Precision amp2117(maxamps), jamp2117(0:maxamps)
      common/to_Ramps_117/amp2117,jamp2117
      Double Precision amp2118(maxamps), jamp2118(0:maxamps)
      common/to_Ramps_118/amp2118,jamp2118
      Double Precision amp2119(maxamps), jamp2119(0:maxamps)
      common/to_Ramps_119/amp2119,jamp2119
      Double Precision amp2120(maxamps), jamp2120(0:maxamps)
      common/to_Ramps_120/amp2120,jamp2120
      Double Precision amp2121(maxamps), jamp2121(0:maxamps)
      common/to_Ramps_121/amp2121,jamp2121
      Double Precision amp2122(maxamps), jamp2122(0:maxamps)
      common/to_Ramps_122/amp2122,jamp2122
      Double Precision amp2123(maxamps), jamp2123(0:maxamps)
      common/to_Ramps_123/amp2123,jamp2123
      Double Precision amp2124(maxamps), jamp2124(0:maxamps)
      common/to_Ramps_124/amp2124,jamp2124
      Double Precision amp2125(maxamps), jamp2125(0:maxamps)
      common/to_Ramps_125/amp2125,jamp2125
      Double Precision amp2126(maxamps), jamp2126(0:maxamps)
      common/to_Ramps_126/amp2126,jamp2126
      Double Precision amp2127(maxamps), jamp2127(0:maxamps)
      common/to_Ramps_127/amp2127,jamp2127
      Double Precision amp2128(maxamps), jamp2128(0:maxamps)
      common/to_Ramps_128/amp2128,jamp2128
      Double Precision amp2129(maxamps), jamp2129(0:maxamps)
      common/to_Ramps_129/amp2129,jamp2129
      Double Precision amp2130(maxamps), jamp2130(0:maxamps)
      common/to_Ramps_130/amp2130,jamp2130
      Double Precision amp2131(maxamps), jamp2131(0:maxamps)
      common/to_Ramps_131/amp2131,jamp2131
      Double Precision amp2132(maxamps), jamp2132(0:maxamps)
      common/to_Ramps_132/amp2132,jamp2132
      Double Precision amp2133(maxamps), jamp2133(0:maxamps)
      common/to_Ramps_133/amp2133,jamp2133
      Double Precision amp2134(maxamps), jamp2134(0:maxamps)
      common/to_Ramps_134/amp2134,jamp2134
      Double Precision amp2135(maxamps), jamp2135(0:maxamps)
      common/to_Ramps_135/amp2135,jamp2135
      Double Precision amp2136(maxamps), jamp2136(0:maxamps)
      common/to_Ramps_136/amp2136,jamp2136
      Double Precision amp2137(maxamps), jamp2137(0:maxamps)
      common/to_Ramps_137/amp2137,jamp2137
      Double Precision amp2138(maxamps), jamp2138(0:maxamps)
      common/to_Ramps_138/amp2138,jamp2138
      Double Precision amp2139(maxamps), jamp2139(0:maxamps)
      common/to_Ramps_139/amp2139,jamp2139
      Double Precision amp2140(maxamps), jamp2140(0:maxamps)
      common/to_Ramps_140/amp2140,jamp2140
      Double Precision amp2141(maxamps), jamp2141(0:maxamps)
      common/to_Ramps_141/amp2141,jamp2141
      Double Precision amp2142(maxamps), jamp2142(0:maxamps)
      common/to_Ramps_142/amp2142,jamp2142
      Double Precision amp2143(maxamps), jamp2143(0:maxamps)
      common/to_Ramps_143/amp2143,jamp2143
      Double Precision amp2144(maxamps), jamp2144(0:maxamps)
      common/to_Ramps_144/amp2144,jamp2144
      Double Precision amp2145(maxamps), jamp2145(0:maxamps)
      common/to_Ramps_145/amp2145,jamp2145
      Double Precision amp2146(maxamps), jamp2146(0:maxamps)
      common/to_Ramps_146/amp2146,jamp2146
      Double Precision amp2147(maxamps), jamp2147(0:maxamps)
      common/to_Ramps_147/amp2147,jamp2147
      Double Precision amp2148(maxamps), jamp2148(0:maxamps)
      common/to_Ramps_148/amp2148,jamp2148
      Double Precision amp2149(maxamps), jamp2149(0:maxamps)
      common/to_Ramps_149/amp2149,jamp2149
      Double Precision amp2150(maxamps), jamp2150(0:maxamps)
      common/to_Ramps_150/amp2150,jamp2150
      Double Precision amp2151(maxamps), jamp2151(0:maxamps)
      common/to_Ramps_151/amp2151,jamp2151
      Double Precision amp2152(maxamps), jamp2152(0:maxamps)
      common/to_Ramps_152/amp2152,jamp2152
      Double Precision amp2153(maxamps), jamp2153(0:maxamps)
      common/to_Ramps_153/amp2153,jamp2153
      Double Precision amp2154(maxamps), jamp2154(0:maxamps)
      common/to_Ramps_154/amp2154,jamp2154
      Double Precision amp2155(maxamps), jamp2155(0:maxamps)
      common/to_Ramps_155/amp2155,jamp2155
      Double Precision amp2156(maxamps), jamp2156(0:maxamps)
      common/to_Ramps_156/amp2156,jamp2156
      Double Precision amp2157(maxamps), jamp2157(0:maxamps)
      common/to_Ramps_157/amp2157,jamp2157
      Double Precision amp2158(maxamps), jamp2158(0:maxamps)
      common/to_Ramps_158/amp2158,jamp2158
      Double Precision amp2159(maxamps), jamp2159(0:maxamps)
      common/to_Ramps_159/amp2159,jamp2159
      Double Precision amp2160(maxamps), jamp2160(0:maxamps)
      common/to_Ramps_160/amp2160,jamp2160
      Double Precision amp2161(maxamps), jamp2161(0:maxamps)
      common/to_Ramps_161/amp2161,jamp2161
      Double Precision amp2162(maxamps), jamp2162(0:maxamps)
      common/to_Ramps_162/amp2162,jamp2162
      Double Precision amp2163(maxamps), jamp2163(0:maxamps)
      common/to_Ramps_163/amp2163,jamp2163
      Double Precision amp2164(maxamps), jamp2164(0:maxamps)
      common/to_Ramps_164/amp2164,jamp2164
      Double Precision amp2165(maxamps), jamp2165(0:maxamps)
      common/to_Ramps_165/amp2165,jamp2165
      Double Precision amp2166(maxamps), jamp2166(0:maxamps)
      common/to_Ramps_166/amp2166,jamp2166
      Double Precision amp2167(maxamps), jamp2167(0:maxamps)
      common/to_Ramps_167/amp2167,jamp2167
      Double Precision amp2168(maxamps), jamp2168(0:maxamps)
      common/to_Ramps_168/amp2168,jamp2168
      Double Precision amp2169(maxamps), jamp2169(0:maxamps)
      common/to_Ramps_169/amp2169,jamp2169
      Double Precision amp2170(maxamps), jamp2170(0:maxamps)
      common/to_Ramps_170/amp2170,jamp2170
      Double Precision amp2171(maxamps), jamp2171(0:maxamps)
      common/to_Ramps_171/amp2171,jamp2171
      Double Precision amp2172(maxamps), jamp2172(0:maxamps)
      common/to_Ramps_172/amp2172,jamp2172
      Double Precision amp2173(maxamps), jamp2173(0:maxamps)
      common/to_Ramps_173/amp2173,jamp2173
      Double Precision amp2174(maxamps), jamp2174(0:maxamps)
      common/to_Ramps_174/amp2174,jamp2174
      Double Precision amp2175(maxamps), jamp2175(0:maxamps)
      common/to_Ramps_175/amp2175,jamp2175
      Double Precision amp2176(maxamps), jamp2176(0:maxamps)
      common/to_Ramps_176/amp2176,jamp2176
      Double Precision amp2177(maxamps), jamp2177(0:maxamps)
      common/to_Ramps_177/amp2177,jamp2177
      Double Precision amp2178(maxamps), jamp2178(0:maxamps)
      common/to_Ramps_178/amp2178,jamp2178
      Double Precision amp2179(maxamps), jamp2179(0:maxamps)
      common/to_Ramps_179/amp2179,jamp2179
      Double Precision amp2180(maxamps), jamp2180(0:maxamps)
      common/to_Ramps_180/amp2180,jamp2180
      Double Precision amp2181(maxamps), jamp2181(0:maxamps)
      common/to_Ramps_181/amp2181,jamp2181
      Double Precision amp2182(maxamps), jamp2182(0:maxamps)
      common/to_Ramps_182/amp2182,jamp2182
      Double Precision amp2183(maxamps), jamp2183(0:maxamps)
      common/to_Ramps_183/amp2183,jamp2183
      Double Precision amp2184(maxamps), jamp2184(0:maxamps)
      common/to_Ramps_184/amp2184,jamp2184
      Double Precision amp2185(maxamps), jamp2185(0:maxamps)
      common/to_Ramps_185/amp2185,jamp2185
      Double Precision amp2186(maxamps), jamp2186(0:maxamps)
      common/to_Ramps_186/amp2186,jamp2186
      Double Precision amp2187(maxamps), jamp2187(0:maxamps)
      common/to_Ramps_187/amp2187,jamp2187
      Double Precision amp2188(maxamps), jamp2188(0:maxamps)
      common/to_Ramps_188/amp2188,jamp2188
      Double Precision amp2189(maxamps), jamp2189(0:maxamps)
      common/to_Ramps_189/amp2189,jamp2189
      Double Precision amp2190(maxamps), jamp2190(0:maxamps)
      common/to_Ramps_190/amp2190,jamp2190
      Double Precision amp2191(maxamps), jamp2191(0:maxamps)
      common/to_Ramps_191/amp2191,jamp2191
      Double Precision amp2192(maxamps), jamp2192(0:maxamps)
      common/to_Ramps_192/amp2192,jamp2192
      Double Precision amp2193(maxamps), jamp2193(0:maxamps)
      common/to_Ramps_193/amp2193,jamp2193
      Double Precision amp2194(maxamps), jamp2194(0:maxamps)
      common/to_Ramps_194/amp2194,jamp2194
      Double Precision amp2195(maxamps), jamp2195(0:maxamps)
      common/to_Ramps_195/amp2195,jamp2195
      Double Precision amp2196(maxamps), jamp2196(0:maxamps)
      common/to_Ramps_196/amp2196,jamp2196
      Double Precision amp2197(maxamps), jamp2197(0:maxamps)
      common/to_Ramps_197/amp2197,jamp2197
      Double Precision amp2198(maxamps), jamp2198(0:maxamps)
      common/to_Ramps_198/amp2198,jamp2198
      Double Precision amp2199(maxamps), jamp2199(0:maxamps)
      common/to_Ramps_199/amp2199,jamp2199
      Double Precision amp2200(maxamps), jamp2200(0:maxamps)
      common/to_Ramps_200/amp2200,jamp2200
      Double Precision amp2201(maxamps), jamp2201(0:maxamps)
      common/to_Ramps_201/amp2201,jamp2201
      Double Precision amp2202(maxamps), jamp2202(0:maxamps)
      common/to_Ramps_202/amp2202,jamp2202
      Double Precision amp2203(maxamps), jamp2203(0:maxamps)
      common/to_Ramps_203/amp2203,jamp2203
      Double Precision amp2204(maxamps), jamp2204(0:maxamps)
      common/to_Ramps_204/amp2204,jamp2204
      Double Precision amp2205(maxamps), jamp2205(0:maxamps)
      common/to_Ramps_205/amp2205,jamp2205
      Double Precision amp2206(maxamps), jamp2206(0:maxamps)
      common/to_Ramps_206/amp2206,jamp2206
      Double Precision amp2207(maxamps), jamp2207(0:maxamps)
      common/to_Ramps_207/amp2207,jamp2207
      Double Precision amp2208(maxamps), jamp2208(0:maxamps)
      common/to_Ramps_208/amp2208,jamp2208
      Double Precision amp2209(maxamps), jamp2209(0:maxamps)
      common/to_Ramps_209/amp2209,jamp2209
      Double Precision amp2210(maxamps), jamp2210(0:maxamps)
      common/to_Ramps_210/amp2210,jamp2210
      Double Precision amp2211(maxamps), jamp2211(0:maxamps)
      common/to_Ramps_211/amp2211,jamp2211
      Double Precision amp2212(maxamps), jamp2212(0:maxamps)
      common/to_Ramps_212/amp2212,jamp2212
      Double Precision amp2213(maxamps), jamp2213(0:maxamps)
      common/to_Ramps_213/amp2213,jamp2213
      Double Precision amp2214(maxamps), jamp2214(0:maxamps)
      common/to_Ramps_214/amp2214,jamp2214
      Double Precision amp2215(maxamps), jamp2215(0:maxamps)
      common/to_Ramps_215/amp2215,jamp2215
      Double Precision amp2216(maxamps), jamp2216(0:maxamps)
      common/to_Ramps_216/amp2216,jamp2216
      Double Precision amp2217(maxamps), jamp2217(0:maxamps)
      common/to_Ramps_217/amp2217,jamp2217
      Double Precision amp2218(maxamps), jamp2218(0:maxamps)
      common/to_Ramps_218/amp2218,jamp2218
      Double Precision amp2219(maxamps), jamp2219(0:maxamps)
      common/to_Ramps_219/amp2219,jamp2219
      Double Precision amp2220(maxamps), jamp2220(0:maxamps)
      common/to_Ramps_220/amp2220,jamp2220
      Double Precision amp2221(maxamps), jamp2221(0:maxamps)
      common/to_Ramps_221/amp2221,jamp2221
      Double Precision amp2222(maxamps), jamp2222(0:maxamps)
      common/to_Ramps_222/amp2222,jamp2222
      Double Precision amp2223(maxamps), jamp2223(0:maxamps)
      common/to_Ramps_223/amp2223,jamp2223
      Double Precision amp2224(maxamps), jamp2224(0:maxamps)
      common/to_Ramps_224/amp2224,jamp2224
      Double Precision amp2225(maxamps), jamp2225(0:maxamps)
      common/to_Ramps_225/amp2225,jamp2225
      Double Precision amp2226(maxamps), jamp2226(0:maxamps)
      common/to_Ramps_226/amp2226,jamp2226
      Double Precision amp2227(maxamps), jamp2227(0:maxamps)
      common/to_Ramps_227/amp2227,jamp2227
      Double Precision amp2228(maxamps), jamp2228(0:maxamps)
      common/to_Ramps_228/amp2228,jamp2228
      Double Precision amp2229(maxamps), jamp2229(0:maxamps)
      common/to_Ramps_229/amp2229,jamp2229
      Double Precision amp2230(maxamps), jamp2230(0:maxamps)
      common/to_Ramps_230/amp2230,jamp2230
      Double Precision amp2231(maxamps), jamp2231(0:maxamps)
      common/to_Ramps_231/amp2231,jamp2231
      Double Precision amp2232(maxamps), jamp2232(0:maxamps)
      common/to_Ramps_232/amp2232,jamp2232
      Double Precision amp2233(maxamps), jamp2233(0:maxamps)
      common/to_Ramps_233/amp2233,jamp2233
      Double Precision amp2234(maxamps), jamp2234(0:maxamps)
      common/to_Ramps_234/amp2234,jamp2234
      Double Precision amp2235(maxamps), jamp2235(0:maxamps)
      common/to_Ramps_235/amp2235,jamp2235
      Double Precision amp2236(maxamps), jamp2236(0:maxamps)
      common/to_Ramps_236/amp2236,jamp2236
      Double Precision amp2237(maxamps), jamp2237(0:maxamps)
      common/to_Ramps_237/amp2237,jamp2237
      Double Precision amp2238(maxamps), jamp2238(0:maxamps)
      common/to_Ramps_238/amp2238,jamp2238
      Double Precision amp2239(maxamps), jamp2239(0:maxamps)
      common/to_Ramps_239/amp2239,jamp2239
      Double Precision amp2240(maxamps), jamp2240(0:maxamps)
      common/to_Ramps_240/amp2240,jamp2240
      Double Precision amp2241(maxamps), jamp2241(0:maxamps)
      common/to_Ramps_241/amp2241,jamp2241
      Double Precision amp2242(maxamps), jamp2242(0:maxamps)
      common/to_Ramps_242/amp2242,jamp2242
      Double Precision amp2243(maxamps), jamp2243(0:maxamps)
      common/to_Ramps_243/amp2243,jamp2243
      Double Precision amp2244(maxamps), jamp2244(0:maxamps)
      common/to_Ramps_244/amp2244,jamp2244
      Double Precision amp2245(maxamps), jamp2245(0:maxamps)
      common/to_Ramps_245/amp2245,jamp2245
      Double Precision amp2246(maxamps), jamp2246(0:maxamps)
      common/to_Ramps_246/amp2246,jamp2246
      Double Precision amp2247(maxamps), jamp2247(0:maxamps)
      common/to_Ramps_247/amp2247,jamp2247
      Double Precision amp2248(maxamps), jamp2248(0:maxamps)
      common/to_Ramps_248/amp2248,jamp2248
      Double Precision amp2249(maxamps), jamp2249(0:maxamps)
      common/to_Ramps_249/amp2249,jamp2249
      Double Precision amp2250(maxamps), jamp2250(0:maxamps)
      common/to_Ramps_250/amp2250,jamp2250
      Double Precision amp2251(maxamps), jamp2251(0:maxamps)
      common/to_Ramps_251/amp2251,jamp2251
      Double Precision amp2252(maxamps), jamp2252(0:maxamps)
      common/to_Ramps_252/amp2252,jamp2252
      Double Precision amp2253(maxamps), jamp2253(0:maxamps)
      common/to_Ramps_253/amp2253,jamp2253
      Double Precision amp2254(maxamps), jamp2254(0:maxamps)
      common/to_Ramps_254/amp2254,jamp2254
      Double Precision amp2255(maxamps), jamp2255(0:maxamps)
      common/to_Ramps_255/amp2255,jamp2255
      Double Precision amp2256(maxamps), jamp2256(0:maxamps)
      common/to_Ramps_256/amp2256,jamp2256
      Double Precision amp2257(maxamps), jamp2257(0:maxamps)
      common/to_Ramps_257/amp2257,jamp2257
      Double Precision amp2258(maxamps), jamp2258(0:maxamps)
      common/to_Ramps_258/amp2258,jamp2258
      Double Precision amp2259(maxamps), jamp2259(0:maxamps)
      common/to_Ramps_259/amp2259,jamp2259
      Double Precision amp2260(maxamps), jamp2260(0:maxamps)
      common/to_Ramps_260/amp2260,jamp2260
      Double Precision amp2261(maxamps), jamp2261(0:maxamps)
      common/to_Ramps_261/amp2261,jamp2261
      Double Precision amp2262(maxamps), jamp2262(0:maxamps)
      common/to_Ramps_262/amp2262,jamp2262
      Double Precision amp2263(maxamps), jamp2263(0:maxamps)
      common/to_Ramps_263/amp2263,jamp2263
      Double Precision amp2264(maxamps), jamp2264(0:maxamps)
      common/to_Ramps_264/amp2264,jamp2264
      Double Precision amp2265(maxamps), jamp2265(0:maxamps)
      common/to_Ramps_265/amp2265,jamp2265
      Double Precision amp2266(maxamps), jamp2266(0:maxamps)
      common/to_Ramps_266/amp2266,jamp2266
      Double Precision amp2267(maxamps), jamp2267(0:maxamps)
      common/to_Ramps_267/amp2267,jamp2267
      Double Precision amp2268(maxamps), jamp2268(0:maxamps)
      common/to_Ramps_268/amp2268,jamp2268
      Double Precision amp2269(maxamps), jamp2269(0:maxamps)
      common/to_Ramps_269/amp2269,jamp2269
      Double Precision amp2270(maxamps), jamp2270(0:maxamps)
      common/to_Ramps_270/amp2270,jamp2270
      Double Precision amp2271(maxamps), jamp2271(0:maxamps)
      common/to_Ramps_271/amp2271,jamp2271
      Double Precision amp2272(maxamps), jamp2272(0:maxamps)
      common/to_Ramps_272/amp2272,jamp2272
      Double Precision amp2273(maxamps), jamp2273(0:maxamps)
      common/to_Ramps_273/amp2273,jamp2273
      Double Precision amp2274(maxamps), jamp2274(0:maxamps)
      common/to_Ramps_274/amp2274,jamp2274
      Double Precision amp2275(maxamps), jamp2275(0:maxamps)
      common/to_Ramps_275/amp2275,jamp2275
      Double Precision amp2276(maxamps), jamp2276(0:maxamps)
      common/to_Ramps_276/amp2276,jamp2276
      Double Precision amp2277(maxamps), jamp2277(0:maxamps)
      common/to_Ramps_277/amp2277,jamp2277
      Double Precision amp2278(maxamps), jamp2278(0:maxamps)
      common/to_Ramps_278/amp2278,jamp2278
      Double Precision amp2279(maxamps), jamp2279(0:maxamps)
      common/to_Ramps_279/amp2279,jamp2279
      Double Precision amp2280(maxamps), jamp2280(0:maxamps)
      common/to_Ramps_280/amp2280,jamp2280
      Double Precision amp2281(maxamps), jamp2281(0:maxamps)
      common/to_Ramps_281/amp2281,jamp2281
      Double Precision amp2282(maxamps), jamp2282(0:maxamps)
      common/to_Ramps_282/amp2282,jamp2282
      Double Precision amp2283(maxamps), jamp2283(0:maxamps)
      common/to_Ramps_283/amp2283,jamp2283
      Double Precision amp2284(maxamps), jamp2284(0:maxamps)
      common/to_Ramps_284/amp2284,jamp2284
      Double Precision amp2285(maxamps), jamp2285(0:maxamps)
      common/to_Ramps_285/amp2285,jamp2285
      Double Precision amp2286(maxamps), jamp2286(0:maxamps)
      common/to_Ramps_286/amp2286,jamp2286
      Double Precision amp2287(maxamps), jamp2287(0:maxamps)
      common/to_Ramps_287/amp2287,jamp2287
      Double Precision amp2288(maxamps), jamp2288(0:maxamps)
      common/to_Ramps_288/amp2288,jamp2288
      Double Precision amp2289(maxamps), jamp2289(0:maxamps)
      common/to_Ramps_289/amp2289,jamp2289
      Double Precision amp2290(maxamps), jamp2290(0:maxamps)
      common/to_Ramps_290/amp2290,jamp2290
      Double Precision amp2291(maxamps), jamp2291(0:maxamps)
      common/to_Ramps_291/amp2291,jamp2291
      Double Precision amp2292(maxamps), jamp2292(0:maxamps)
      common/to_Ramps_292/amp2292,jamp2292
      Double Precision amp2293(maxamps), jamp2293(0:maxamps)
      common/to_Ramps_293/amp2293,jamp2293
      Double Precision amp2294(maxamps), jamp2294(0:maxamps)
      common/to_Ramps_294/amp2294,jamp2294
      Double Precision amp2295(maxamps), jamp2295(0:maxamps)
      common/to_Ramps_295/amp2295,jamp2295
      Double Precision amp2296(maxamps), jamp2296(0:maxamps)
      common/to_Ramps_296/amp2296,jamp2296
      Double Precision amp2297(maxamps), jamp2297(0:maxamps)
      common/to_Ramps_297/amp2297,jamp2297
      Double Precision amp2298(maxamps), jamp2298(0:maxamps)
      common/to_Ramps_298/amp2298,jamp2298
      Double Precision amp2299(maxamps), jamp2299(0:maxamps)
      common/to_Ramps_299/amp2299,jamp2299
      Double Precision amp2300(maxamps), jamp2300(0:maxamps)
      common/to_Ramps_300/amp2300,jamp2300
      Double Precision amp2301(maxamps), jamp2301(0:maxamps)
      common/to_Ramps_301/amp2301,jamp2301
      Double Precision amp2302(maxamps), jamp2302(0:maxamps)
      common/to_Ramps_302/amp2302,jamp2302
      Double Precision amp2303(maxamps), jamp2303(0:maxamps)
      common/to_Ramps_303/amp2303,jamp2303
      Double Precision amp2304(maxamps), jamp2304(0:maxamps)
      common/to_Ramps_304/amp2304,jamp2304
      Double Precision amp2305(maxamps), jamp2305(0:maxamps)
      common/to_Ramps_305/amp2305,jamp2305
      Double Precision amp2306(maxamps), jamp2306(0:maxamps)
      common/to_Ramps_306/amp2306,jamp2306
      Double Precision amp2307(maxamps), jamp2307(0:maxamps)
      common/to_Ramps_307/amp2307,jamp2307
      Double Precision amp2308(maxamps), jamp2308(0:maxamps)
      common/to_Ramps_308/amp2308,jamp2308
      Double Precision amp2309(maxamps), jamp2309(0:maxamps)
      common/to_Ramps_309/amp2309,jamp2309
      Double Precision amp2310(maxamps), jamp2310(0:maxamps)
      common/to_Ramps_310/amp2310,jamp2310
      Double Precision amp2311(maxamps), jamp2311(0:maxamps)
      common/to_Ramps_311/amp2311,jamp2311
      Double Precision amp2312(maxamps), jamp2312(0:maxamps)
      common/to_Ramps_312/amp2312,jamp2312
      Double Precision amp2313(maxamps), jamp2313(0:maxamps)
      common/to_Ramps_313/amp2313,jamp2313
      Double Precision amp2314(maxamps), jamp2314(0:maxamps)
      common/to_Ramps_314/amp2314,jamp2314
      Double Precision amp2315(maxamps), jamp2315(0:maxamps)
      common/to_Ramps_315/amp2315,jamp2315
      Double Precision amp2316(maxamps), jamp2316(0:maxamps)
      common/to_Ramps_316/amp2316,jamp2316
      Double Precision amp2317(maxamps), jamp2317(0:maxamps)
      common/to_Ramps_317/amp2317,jamp2317
      Double Precision amp2318(maxamps), jamp2318(0:maxamps)
      common/to_Ramps_318/amp2318,jamp2318
      Double Precision amp2319(maxamps), jamp2319(0:maxamps)
      common/to_Ramps_319/amp2319,jamp2319
      Double Precision amp2320(maxamps), jamp2320(0:maxamps)
      common/to_Ramps_320/amp2320,jamp2320
      Double Precision amp2321(maxamps), jamp2321(0:maxamps)
      common/to_Ramps_321/amp2321,jamp2321
      Double Precision amp2322(maxamps), jamp2322(0:maxamps)
      common/to_Ramps_322/amp2322,jamp2322
      Double Precision amp2323(maxamps), jamp2323(0:maxamps)
      common/to_Ramps_323/amp2323,jamp2323
      Double Precision amp2324(maxamps), jamp2324(0:maxamps)
      common/to_Ramps_324/amp2324,jamp2324
      Double Precision amp2325(maxamps), jamp2325(0:maxamps)
      common/to_Ramps_325/amp2325,jamp2325
      Double Precision amp2326(maxamps), jamp2326(0:maxamps)
      common/to_Ramps_326/amp2326,jamp2326
      Double Precision amp2327(maxamps), jamp2327(0:maxamps)
      common/to_Ramps_327/amp2327,jamp2327
      Double Precision amp2328(maxamps), jamp2328(0:maxamps)
      common/to_Ramps_328/amp2328,jamp2328
      Double Precision amp2329(maxamps), jamp2329(0:maxamps)
      common/to_Ramps_329/amp2329,jamp2329
      Double Precision amp2330(maxamps), jamp2330(0:maxamps)
      common/to_Ramps_330/amp2330,jamp2330
      Double Precision amp2331(maxamps), jamp2331(0:maxamps)
      common/to_Ramps_331/amp2331,jamp2331
      Double Precision amp2332(maxamps), jamp2332(0:maxamps)
      common/to_Ramps_332/amp2332,jamp2332
      Double Precision amp2333(maxamps), jamp2333(0:maxamps)
      common/to_Ramps_333/amp2333,jamp2333
      Double Precision amp2334(maxamps), jamp2334(0:maxamps)
      common/to_Ramps_334/amp2334,jamp2334
      Double Precision amp2335(maxamps), jamp2335(0:maxamps)
      common/to_Ramps_335/amp2335,jamp2335
      Double Precision amp2336(maxamps), jamp2336(0:maxamps)
      common/to_Ramps_336/amp2336,jamp2336
      Double Precision amp2337(maxamps), jamp2337(0:maxamps)
      common/to_Ramps_337/amp2337,jamp2337
      Double Precision amp2338(maxamps), jamp2338(0:maxamps)
      common/to_Ramps_338/amp2338,jamp2338
      Double Precision amp2339(maxamps), jamp2339(0:maxamps)
      common/to_Ramps_339/amp2339,jamp2339
      Double Precision amp2340(maxamps), jamp2340(0:maxamps)
      common/to_Ramps_340/amp2340,jamp2340
      Double Precision amp2341(maxamps), jamp2341(0:maxamps)
      common/to_Ramps_341/amp2341,jamp2341
      Double Precision amp2342(maxamps), jamp2342(0:maxamps)
      common/to_Ramps_342/amp2342,jamp2342
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal,maxamps)
      integer color(2,nexternal),color1(2,nexternal)
      double precision random,xtarget
      external random
      integer legs(nexternal),lstr,i,j
      character*140 str
      integer ic(nexternal),legs1(nexternal)
      integer iflow,ifl
      logical mtc,even
      
      do i=1,nexternal
         ic(i)=i
      enddo
      mtc=.false.
 10   call nexper(nexternal- 5,ic( 5+1),mtc,even)
      do i= 5+1,nexternal
         ic(i)=ic(i)+ 5
      enddo
      CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL)
      
      call convert_to_string(nexternal,legs1,str,lstr)
      
      if (str.eq."-1-125-1112-1-2") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif (str.eq."-1-125-1112-1-4") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif (str.eq."-1125-11121-2") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif (str.eq."-1125-11121-4") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif (str.eq."-1125-1112-23") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif (str.eq."-1125-1112-25") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif (str.eq."-1125-1112-43") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif (str.eq."-1125-1112-45") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif (str.eq."-1-225-1112-2-2") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif (str.eq."-1-225-1112-2-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif (str.eq."-1225-11121-1") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif (str.eq."-1225-1112-13") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif (str.eq."-1225-1112-15") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif (str.eq."-1225-11122-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif (str.eq."-1225-11122-4") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif (str.eq."-1225-11124-4") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif (str.eq."-1225-11123-3") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif (str.eq."-1225-11125-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif (str.eq."-1225-111200") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif (str.eq."-1-425-1112-2-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif (str.eq."-1-425-1112-4-4") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif (str.eq."-1425-11121-1") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif (str.eq."-1425-1112-13") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif (str.eq."-1425-1112-15") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif (str.eq."-1425-11122-2") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif (str.eq."-1425-1112-24") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif (str.eq."-1425-11124-4") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif (str.eq."-1425-11123-3") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif (str.eq."-1425-11125-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif (str.eq."-1425-111200") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
      elseif (str.eq."-1-325-1112-1-2") then
         include "leshouches_R_031.inc"
         iflow=nint(jamp2031(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2031(i)
         enddo
         goto 20
      elseif (str.eq."-1-325-1112-1-4") then
         include "leshouches_R_032.inc"
         iflow=nint(jamp2032(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2032(i)
         enddo
         goto 20
      elseif (str.eq."-1-325-1112-2-3") then
         include "leshouches_R_033.inc"
         iflow=nint(jamp2033(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2033(i)
         enddo
         goto 20
      elseif (str.eq."-1-325-1112-4-3") then
         include "leshouches_R_034.inc"
         iflow=nint(jamp2034(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2034(i)
         enddo
         goto 20
      elseif (str.eq."-1325-1112-23") then
         include "leshouches_R_035.inc"
         iflow=nint(jamp2035(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2035(i)
         enddo
         goto 20
      elseif (str.eq."-1325-1112-43") then
         include "leshouches_R_036.inc"
         iflow=nint(jamp2036(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2036(i)
         enddo
         goto 20
      elseif (str.eq."-1-525-1112-1-2") then
         include "leshouches_R_037.inc"
         iflow=nint(jamp2037(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2037(i)
         enddo
         goto 20
      elseif (str.eq."-1-525-1112-1-4") then
         include "leshouches_R_038.inc"
         iflow=nint(jamp2038(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2038(i)
         enddo
         goto 20
      elseif (str.eq."-1-525-1112-2-5") then
         include "leshouches_R_039.inc"
         iflow=nint(jamp2039(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2039(i)
         enddo
         goto 20
      elseif (str.eq."-1-525-1112-4-5") then
         include "leshouches_R_040.inc"
         iflow=nint(jamp2040(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2040(i)
         enddo
         goto 20
      elseif (str.eq."-1525-1112-25") then
         include "leshouches_R_041.inc"
         iflow=nint(jamp2041(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2041(i)
         enddo
         goto 20
      elseif (str.eq."-1525-1112-45") then
         include "leshouches_R_042.inc"
         iflow=nint(jamp2042(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2042(i)
         enddo
         goto 20
      elseif (str.eq."-1025-1112-20") then
         include "leshouches_R_043.inc"
         iflow=nint(jamp2043(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2043(i)
         enddo
         goto 20
      elseif (str.eq."-1025-1112-40") then
         include "leshouches_R_044.inc"
         iflow=nint(jamp2044(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2044(i)
         enddo
         goto 20
      elseif (str.eq."1-125-11121-2") then
         include "leshouches_R_045.inc"
         iflow=nint(jamp2045(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2045(i)
         enddo
         goto 20
      elseif (str.eq."1-125-11121-4") then
         include "leshouches_R_046.inc"
         iflow=nint(jamp2046(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2046(i)
         enddo
         goto 20
      elseif (str.eq."1-125-1112-23") then
         include "leshouches_R_047.inc"
         iflow=nint(jamp2047(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2047(i)
         enddo
         goto 20
      elseif (str.eq."1-125-1112-25") then
         include "leshouches_R_048.inc"
         iflow=nint(jamp2048(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2048(i)
         enddo
         goto 20
      elseif (str.eq."1-125-1112-43") then
         include "leshouches_R_049.inc"
         iflow=nint(jamp2049(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2049(i)
         enddo
         goto 20
      elseif (str.eq."1-125-1112-45") then
         include "leshouches_R_050.inc"
         iflow=nint(jamp2050(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2050(i)
         enddo
         goto 20
      elseif (str.eq."1225-111211") then
         include "leshouches_R_051.inc"
         iflow=nint(jamp2051(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2051(i)
         enddo
         goto 20
      elseif (str.eq."1225-111213") then
         include "leshouches_R_052.inc"
         iflow=nint(jamp2052(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2052(i)
         enddo
         goto 20
      elseif (str.eq."1225-111215") then
         include "leshouches_R_053.inc"
         iflow=nint(jamp2053(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2053(i)
         enddo
         goto 20
      elseif (str.eq."1425-111211") then
         include "leshouches_R_054.inc"
         iflow=nint(jamp2054(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2054(i)
         enddo
         goto 20
      elseif (str.eq."1425-111213") then
         include "leshouches_R_055.inc"
         iflow=nint(jamp2055(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2055(i)
         enddo
         goto 20
      elseif (str.eq."1425-111215") then
         include "leshouches_R_056.inc"
         iflow=nint(jamp2056(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2056(i)
         enddo
         goto 20
      elseif (str.eq."1-325-11121-2") then
         include "leshouches_R_057.inc"
         iflow=nint(jamp2057(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2057(i)
         enddo
         goto 20
      elseif (str.eq."1-325-11121-4") then
         include "leshouches_R_058.inc"
         iflow=nint(jamp2058(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2058(i)
         enddo
         goto 20
      elseif (str.eq."1-525-11121-2") then
         include "leshouches_R_059.inc"
         iflow=nint(jamp2059(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2059(i)
         enddo
         goto 20
      elseif (str.eq."1-525-11121-4") then
         include "leshouches_R_060.inc"
         iflow=nint(jamp2060(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2060(i)
         enddo
         goto 20
      elseif (str.eq."-2-125-1112-2-2") then
         include "leshouches_R_061.inc"
         iflow=nint(jamp2061(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2061(i)
         enddo
         goto 20
      elseif (str.eq."-2-125-1112-2-4") then
         include "leshouches_R_062.inc"
         iflow=nint(jamp2062(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2062(i)
         enddo
         goto 20
      elseif (str.eq."-2225-11121-2") then
         include "leshouches_R_063.inc"
         iflow=nint(jamp2063(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2063(i)
         enddo
         goto 20
      elseif (str.eq."-2225-11121-4") then
         include "leshouches_R_064.inc"
         iflow=nint(jamp2064(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2064(i)
         enddo
         goto 20
      elseif (str.eq."-2225-1112-23") then
         include "leshouches_R_065.inc"
         iflow=nint(jamp2065(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2065(i)
         enddo
         goto 20
      elseif (str.eq."-2225-1112-25") then
         include "leshouches_R_066.inc"
         iflow=nint(jamp2066(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2066(i)
         enddo
         goto 20
      elseif (str.eq."-2225-1112-43") then
         include "leshouches_R_067.inc"
         iflow=nint(jamp2067(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2067(i)
         enddo
         goto 20
      elseif (str.eq."-2225-1112-45") then
         include "leshouches_R_068.inc"
         iflow=nint(jamp2068(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2068(i)
         enddo
         goto 20
      elseif (str.eq."-2425-11121-2") then
         include "leshouches_R_069.inc"
         iflow=nint(jamp2069(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2069(i)
         enddo
         goto 20
      elseif (str.eq."-2425-1112-23") then
         include "leshouches_R_070.inc"
         iflow=nint(jamp2070(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2070(i)
         enddo
         goto 20
      elseif (str.eq."-2425-1112-25") then
         include "leshouches_R_071.inc"
         iflow=nint(jamp2071(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2071(i)
         enddo
         goto 20
      elseif (str.eq."-2-325-1112-2-2") then
         include "leshouches_R_072.inc"
         iflow=nint(jamp2072(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2072(i)
         enddo
         goto 20
      elseif (str.eq."-2-325-1112-2-4") then
         include "leshouches_R_073.inc"
         iflow=nint(jamp2073(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2073(i)
         enddo
         goto 20
      elseif (str.eq."-2-525-1112-2-2") then
         include "leshouches_R_074.inc"
         iflow=nint(jamp2074(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2074(i)
         enddo
         goto 20
      elseif (str.eq."-2-525-1112-2-4") then
         include "leshouches_R_075.inc"
         iflow=nint(jamp2075(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2075(i)
         enddo
         goto 20
      elseif (str.eq."2-125-11121-1") then
         include "leshouches_R_076.inc"
         iflow=nint(jamp2076(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2076(i)
         enddo
         goto 20
      elseif (str.eq."2-125-1112-13") then
         include "leshouches_R_077.inc"
         iflow=nint(jamp2077(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2077(i)
         enddo
         goto 20
      elseif (str.eq."2-125-1112-15") then
         include "leshouches_R_078.inc"
         iflow=nint(jamp2078(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2078(i)
         enddo
         goto 20
      elseif (str.eq."2-125-11122-2") then
         include "leshouches_R_079.inc"
         iflow=nint(jamp2079(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2079(i)
         enddo
         goto 20
      elseif (str.eq."2-125-11122-4") then
         include "leshouches_R_080.inc"
         iflow=nint(jamp2080(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2080(i)
         enddo
         goto 20
      elseif (str.eq."2-125-11124-4") then
         include "leshouches_R_081.inc"
         iflow=nint(jamp2081(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2081(i)
         enddo
         goto 20
      elseif (str.eq."2-125-11123-3") then
         include "leshouches_R_082.inc"
         iflow=nint(jamp2082(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2082(i)
         enddo
         goto 20
      elseif (str.eq."2-125-11125-5") then
         include "leshouches_R_083.inc"
         iflow=nint(jamp2083(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2083(i)
         enddo
         goto 20
      elseif (str.eq."2-125-111200") then
         include "leshouches_R_084.inc"
         iflow=nint(jamp2084(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2084(i)
         enddo
         goto 20
      elseif (str.eq."2125-111211") then
         include "leshouches_R_085.inc"
         iflow=nint(jamp2085(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2085(i)
         enddo
         goto 20
      elseif (str.eq."2125-111213") then
         include "leshouches_R_086.inc"
         iflow=nint(jamp2086(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2086(i)
         enddo
         goto 20
      elseif (str.eq."2125-111215") then
         include "leshouches_R_087.inc"
         iflow=nint(jamp2087(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2087(i)
         enddo
         goto 20
      elseif (str.eq."2-225-11121-2") then
         include "leshouches_R_088.inc"
         iflow=nint(jamp2088(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2088(i)
         enddo
         goto 20
      elseif (str.eq."2-225-11121-4") then
         include "leshouches_R_089.inc"
         iflow=nint(jamp2089(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2089(i)
         enddo
         goto 20
      elseif (str.eq."2-225-1112-23") then
         include "leshouches_R_090.inc"
         iflow=nint(jamp2090(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2090(i)
         enddo
         goto 20
      elseif (str.eq."2-225-1112-25") then
         include "leshouches_R_091.inc"
         iflow=nint(jamp2091(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2091(i)
         enddo
         goto 20
      elseif (str.eq."2-225-1112-43") then
         include "leshouches_R_092.inc"
         iflow=nint(jamp2092(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2092(i)
         enddo
         goto 20
      elseif (str.eq."2-225-1112-45") then
         include "leshouches_R_093.inc"
         iflow=nint(jamp2093(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2093(i)
         enddo
         goto 20
      elseif (str.eq."2225-111212") then
         include "leshouches_R_094.inc"
         iflow=nint(jamp2094(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2094(i)
         enddo
         goto 20
      elseif (str.eq."2225-111223") then
         include "leshouches_R_095.inc"
         iflow=nint(jamp2095(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2095(i)
         enddo
         goto 20
      elseif (str.eq."2225-111225") then
         include "leshouches_R_096.inc"
         iflow=nint(jamp2096(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2096(i)
         enddo
         goto 20
      elseif (str.eq."2-425-11121-4") then
         include "leshouches_R_097.inc"
         iflow=nint(jamp2097(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2097(i)
         enddo
         goto 20
      elseif (str.eq."2-425-1112-43") then
         include "leshouches_R_098.inc"
         iflow=nint(jamp2098(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2098(i)
         enddo
         goto 20
      elseif (str.eq."2-425-1112-45") then
         include "leshouches_R_099.inc"
         iflow=nint(jamp2099(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2099(i)
         enddo
         goto 20
      elseif (str.eq."2425-111212") then
         include "leshouches_R_100.inc"
         iflow=nint(jamp2100(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2100(i)
         enddo
         goto 20
      elseif (str.eq."2425-111214") then
         include "leshouches_R_101.inc"
         iflow=nint(jamp2101(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2101(i)
         enddo
         goto 20
      elseif (str.eq."2425-111223") then
         include "leshouches_R_102.inc"
         iflow=nint(jamp2102(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2102(i)
         enddo
         goto 20
      elseif (str.eq."2425-111225") then
         include "leshouches_R_103.inc"
         iflow=nint(jamp2103(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2103(i)
         enddo
         goto 20
      elseif (str.eq."2425-111243") then
         include "leshouches_R_104.inc"
         iflow=nint(jamp2104(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2104(i)
         enddo
         goto 20
      elseif (str.eq."2425-111245") then
         include "leshouches_R_105.inc"
         iflow=nint(jamp2105(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2105(i)
         enddo
         goto 20
      elseif (str.eq."2-325-11121-1") then
         include "leshouches_R_106.inc"
         iflow=nint(jamp2106(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2106(i)
         enddo
         goto 20
      elseif (str.eq."2-325-11121-3") then
         include "leshouches_R_107.inc"
         iflow=nint(jamp2107(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2107(i)
         enddo
         goto 20
      elseif (str.eq."2-325-11122-2") then
         include "leshouches_R_108.inc"
         iflow=nint(jamp2108(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2108(i)
         enddo
         goto 20
      elseif (str.eq."2-325-11122-4") then
         include "leshouches_R_109.inc"
         iflow=nint(jamp2109(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2109(i)
         enddo
         goto 20
      elseif (str.eq."2-325-11124-4") then
         include "leshouches_R_110.inc"
         iflow=nint(jamp2110(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2110(i)
         enddo
         goto 20
      elseif (str.eq."2-325-11123-3") then
         include "leshouches_R_111.inc"
         iflow=nint(jamp2111(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2111(i)
         enddo
         goto 20
      elseif (str.eq."2-325-1112-35") then
         include "leshouches_R_112.inc"
         iflow=nint(jamp2112(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2112(i)
         enddo
         goto 20
      elseif (str.eq."2-325-11125-5") then
         include "leshouches_R_113.inc"
         iflow=nint(jamp2113(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2113(i)
         enddo
         goto 20
      elseif (str.eq."2-325-111200") then
         include "leshouches_R_114.inc"
         iflow=nint(jamp2114(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2114(i)
         enddo
         goto 20
      elseif (str.eq."2325-111213") then
         include "leshouches_R_115.inc"
         iflow=nint(jamp2115(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2115(i)
         enddo
         goto 20
      elseif (str.eq."2325-111233") then
         include "leshouches_R_116.inc"
         iflow=nint(jamp2116(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2116(i)
         enddo
         goto 20
      elseif (str.eq."2325-111235") then
         include "leshouches_R_117.inc"
         iflow=nint(jamp2117(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2117(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11121-1") then
         include "leshouches_R_118.inc"
         iflow=nint(jamp2118(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2118(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11121-5") then
         include "leshouches_R_119.inc"
         iflow=nint(jamp2119(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2119(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11122-2") then
         include "leshouches_R_120.inc"
         iflow=nint(jamp2120(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2120(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11122-4") then
         include "leshouches_R_121.inc"
         iflow=nint(jamp2121(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2121(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11124-4") then
         include "leshouches_R_122.inc"
         iflow=nint(jamp2122(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2122(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11123-3") then
         include "leshouches_R_123.inc"
         iflow=nint(jamp2123(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2123(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11123-5") then
         include "leshouches_R_124.inc"
         iflow=nint(jamp2124(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2124(i)
         enddo
         goto 20
      elseif (str.eq."2-525-11125-5") then
         include "leshouches_R_125.inc"
         iflow=nint(jamp2125(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2125(i)
         enddo
         goto 20
      elseif (str.eq."2-525-111200") then
         include "leshouches_R_126.inc"
         iflow=nint(jamp2126(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2126(i)
         enddo
         goto 20
      elseif (str.eq."2525-111215") then
         include "leshouches_R_127.inc"
         iflow=nint(jamp2127(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2127(i)
         enddo
         goto 20
      elseif (str.eq."2525-111235") then
         include "leshouches_R_128.inc"
         iflow=nint(jamp2128(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2128(i)
         enddo
         goto 20
      elseif (str.eq."2525-111255") then
         include "leshouches_R_129.inc"
         iflow=nint(jamp2129(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2129(i)
         enddo
         goto 20
      elseif (str.eq."2025-111210") then
         include "leshouches_R_130.inc"
         iflow=nint(jamp2130(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2130(i)
         enddo
         goto 20
      elseif (str.eq."2025-111230") then
         include "leshouches_R_131.inc"
         iflow=nint(jamp2131(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2131(i)
         enddo
         goto 20
      elseif (str.eq."2025-111250") then
         include "leshouches_R_132.inc"
         iflow=nint(jamp2132(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2132(i)
         enddo
         goto 20
      elseif (str.eq."-4-125-1112-2-4") then
         include "leshouches_R_133.inc"
         iflow=nint(jamp2133(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2133(i)
         enddo
         goto 20
      elseif (str.eq."-4-125-1112-4-4") then
         include "leshouches_R_134.inc"
         iflow=nint(jamp2134(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2134(i)
         enddo
         goto 20
      elseif (str.eq."-4225-11121-4") then
         include "leshouches_R_135.inc"
         iflow=nint(jamp2135(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2135(i)
         enddo
         goto 20
      elseif (str.eq."-4225-1112-43") then
         include "leshouches_R_136.inc"
         iflow=nint(jamp2136(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2136(i)
         enddo
         goto 20
      elseif (str.eq."-4225-1112-45") then
         include "leshouches_R_137.inc"
         iflow=nint(jamp2137(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2137(i)
         enddo
         goto 20
      elseif (str.eq."-4425-11121-2") then
         include "leshouches_R_138.inc"
         iflow=nint(jamp2138(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2138(i)
         enddo
         goto 20
      elseif (str.eq."-4425-11121-4") then
         include "leshouches_R_139.inc"
         iflow=nint(jamp2139(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2139(i)
         enddo
         goto 20
      elseif (str.eq."-4425-1112-23") then
         include "leshouches_R_140.inc"
         iflow=nint(jamp2140(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2140(i)
         enddo
         goto 20
      elseif (str.eq."-4425-1112-25") then
         include "leshouches_R_141.inc"
         iflow=nint(jamp2141(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2141(i)
         enddo
         goto 20
      elseif (str.eq."-4425-1112-43") then
         include "leshouches_R_142.inc"
         iflow=nint(jamp2142(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2142(i)
         enddo
         goto 20
      elseif (str.eq."-4425-1112-45") then
         include "leshouches_R_143.inc"
         iflow=nint(jamp2143(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2143(i)
         enddo
         goto 20
      elseif (str.eq."-4-325-1112-2-4") then
         include "leshouches_R_144.inc"
         iflow=nint(jamp2144(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2144(i)
         enddo
         goto 20
      elseif (str.eq."-4-325-1112-4-4") then
         include "leshouches_R_145.inc"
         iflow=nint(jamp2145(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2145(i)
         enddo
         goto 20
      elseif (str.eq."-4-525-1112-2-4") then
         include "leshouches_R_146.inc"
         iflow=nint(jamp2146(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2146(i)
         enddo
         goto 20
      elseif (str.eq."-4-525-1112-4-4") then
         include "leshouches_R_147.inc"
         iflow=nint(jamp2147(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2147(i)
         enddo
         goto 20
      elseif (str.eq."4-125-11121-1") then
         include "leshouches_R_148.inc"
         iflow=nint(jamp2148(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2148(i)
         enddo
         goto 20
      elseif (str.eq."4-125-1112-13") then
         include "leshouches_R_149.inc"
         iflow=nint(jamp2149(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2149(i)
         enddo
         goto 20
      elseif (str.eq."4-125-1112-15") then
         include "leshouches_R_150.inc"
         iflow=nint(jamp2150(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2150(i)
         enddo
         goto 20
      elseif (str.eq."4-125-11122-2") then
         include "leshouches_R_151.inc"
         iflow=nint(jamp2151(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2151(i)
         enddo
         goto 20
      elseif (str.eq."4-125-1112-24") then
         include "leshouches_R_152.inc"
         iflow=nint(jamp2152(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2152(i)
         enddo
         goto 20
      elseif (str.eq."4-125-11124-4") then
         include "leshouches_R_153.inc"
         iflow=nint(jamp2153(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2153(i)
         enddo
         goto 20
      elseif (str.eq."4-125-11123-3") then
         include "leshouches_R_154.inc"
         iflow=nint(jamp2154(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2154(i)
         enddo
         goto 20
      elseif (str.eq."4-125-11125-5") then
         include "leshouches_R_155.inc"
         iflow=nint(jamp2155(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2155(i)
         enddo
         goto 20
      elseif (str.eq."4-125-111200") then
         include "leshouches_R_156.inc"
         iflow=nint(jamp2156(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2156(i)
         enddo
         goto 20
      elseif (str.eq."4125-111211") then
         include "leshouches_R_157.inc"
         iflow=nint(jamp2157(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2157(i)
         enddo
         goto 20
      elseif (str.eq."4125-111213") then
         include "leshouches_R_158.inc"
         iflow=nint(jamp2158(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2158(i)
         enddo
         goto 20
      elseif (str.eq."4125-111215") then
         include "leshouches_R_159.inc"
         iflow=nint(jamp2159(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2159(i)
         enddo
         goto 20
      elseif (str.eq."4-225-11121-2") then
         include "leshouches_R_160.inc"
         iflow=nint(jamp2160(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2160(i)
         enddo
         goto 20
      elseif (str.eq."4-225-1112-23") then
         include "leshouches_R_161.inc"
         iflow=nint(jamp2161(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2161(i)
         enddo
         goto 20
      elseif (str.eq."4-225-1112-25") then
         include "leshouches_R_162.inc"
         iflow=nint(jamp2162(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2162(i)
         enddo
         goto 20
      elseif (str.eq."4225-111212") then
         include "leshouches_R_163.inc"
         iflow=nint(jamp2163(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2163(i)
         enddo
         goto 20
      elseif (str.eq."4225-111214") then
         include "leshouches_R_164.inc"
         iflow=nint(jamp2164(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2164(i)
         enddo
         goto 20
      elseif (str.eq."4225-111223") then
         include "leshouches_R_165.inc"
         iflow=nint(jamp2165(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2165(i)
         enddo
         goto 20
      elseif (str.eq."4225-111225") then
         include "leshouches_R_166.inc"
         iflow=nint(jamp2166(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2166(i)
         enddo
         goto 20
      elseif (str.eq."4225-111243") then
         include "leshouches_R_167.inc"
         iflow=nint(jamp2167(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2167(i)
         enddo
         goto 20
      elseif (str.eq."4225-111245") then
         include "leshouches_R_168.inc"
         iflow=nint(jamp2168(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2168(i)
         enddo
         goto 20
      elseif (str.eq."4-425-11121-2") then
         include "leshouches_R_169.inc"
         iflow=nint(jamp2169(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2169(i)
         enddo
         goto 20
      elseif (str.eq."4-425-11121-4") then
         include "leshouches_R_170.inc"
         iflow=nint(jamp2170(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2170(i)
         enddo
         goto 20
      elseif (str.eq."4-425-1112-23") then
         include "leshouches_R_171.inc"
         iflow=nint(jamp2171(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2171(i)
         enddo
         goto 20
      elseif (str.eq."4-425-1112-25") then
         include "leshouches_R_172.inc"
         iflow=nint(jamp2172(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2172(i)
         enddo
         goto 20
      elseif (str.eq."4-425-1112-43") then
         include "leshouches_R_173.inc"
         iflow=nint(jamp2173(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2173(i)
         enddo
         goto 20
      elseif (str.eq."4-425-1112-45") then
         include "leshouches_R_174.inc"
         iflow=nint(jamp2174(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2174(i)
         enddo
         goto 20
      elseif (str.eq."4425-111214") then
         include "leshouches_R_175.inc"
         iflow=nint(jamp2175(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2175(i)
         enddo
         goto 20
      elseif (str.eq."4425-111243") then
         include "leshouches_R_176.inc"
         iflow=nint(jamp2176(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2176(i)
         enddo
         goto 20
      elseif (str.eq."4425-111245") then
         include "leshouches_R_177.inc"
         iflow=nint(jamp2177(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2177(i)
         enddo
         goto 20
      elseif (str.eq."4-325-11121-1") then
         include "leshouches_R_178.inc"
         iflow=nint(jamp2178(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2178(i)
         enddo
         goto 20
      elseif (str.eq."4-325-11121-3") then
         include "leshouches_R_179.inc"
         iflow=nint(jamp2179(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2179(i)
         enddo
         goto 20
      elseif (str.eq."4-325-11122-2") then
         include "leshouches_R_180.inc"
         iflow=nint(jamp2180(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2180(i)
         enddo
         goto 20
      elseif (str.eq."4-325-1112-24") then
         include "leshouches_R_181.inc"
         iflow=nint(jamp2181(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2181(i)
         enddo
         goto 20
      elseif (str.eq."4-325-11124-4") then
         include "leshouches_R_182.inc"
         iflow=nint(jamp2182(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2182(i)
         enddo
         goto 20
      elseif (str.eq."4-325-11123-3") then
         include "leshouches_R_183.inc"
         iflow=nint(jamp2183(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2183(i)
         enddo
         goto 20
      elseif (str.eq."4-325-1112-35") then
         include "leshouches_R_184.inc"
         iflow=nint(jamp2184(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2184(i)
         enddo
         goto 20
      elseif (str.eq."4-325-11125-5") then
         include "leshouches_R_185.inc"
         iflow=nint(jamp2185(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2185(i)
         enddo
         goto 20
      elseif (str.eq."4-325-111200") then
         include "leshouches_R_186.inc"
         iflow=nint(jamp2186(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2186(i)
         enddo
         goto 20
      elseif (str.eq."4325-111213") then
         include "leshouches_R_187.inc"
         iflow=nint(jamp2187(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2187(i)
         enddo
         goto 20
      elseif (str.eq."4325-111233") then
         include "leshouches_R_188.inc"
         iflow=nint(jamp2188(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2188(i)
         enddo
         goto 20
      elseif (str.eq."4325-111235") then
         include "leshouches_R_189.inc"
         iflow=nint(jamp2189(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2189(i)
         enddo
         goto 20
      elseif (str.eq."4-525-11121-1") then
         include "leshouches_R_190.inc"
         iflow=nint(jamp2190(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2190(i)
         enddo
         goto 20
      elseif (str.eq."4-525-11121-5") then
         include "leshouches_R_191.inc"
         iflow=nint(jamp2191(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2191(i)
         enddo
         goto 20
      elseif (str.eq."4-525-11122-2") then
         include "leshouches_R_192.inc"
         iflow=nint(jamp2192(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2192(i)
         enddo
         goto 20
      elseif (str.eq."4-525-1112-24") then
         include "leshouches_R_193.inc"
         iflow=nint(jamp2193(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2193(i)
         enddo
         goto 20
      elseif (str.eq."4-525-11124-4") then
         include "leshouches_R_194.inc"
         iflow=nint(jamp2194(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2194(i)
         enddo
         goto 20
      elseif (str.eq."4-525-11123-3") then
         include "leshouches_R_195.inc"
         iflow=nint(jamp2195(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2195(i)
         enddo
         goto 20
      elseif (str.eq."4-525-11123-5") then
         include "leshouches_R_196.inc"
         iflow=nint(jamp2196(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2196(i)
         enddo
         goto 20
      elseif (str.eq."4-525-11125-5") then
         include "leshouches_R_197.inc"
         iflow=nint(jamp2197(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2197(i)
         enddo
         goto 20
      elseif (str.eq."4-525-111200") then
         include "leshouches_R_198.inc"
         iflow=nint(jamp2198(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2198(i)
         enddo
         goto 20
      elseif (str.eq."4525-111215") then
         include "leshouches_R_199.inc"
         iflow=nint(jamp2199(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2199(i)
         enddo
         goto 20
      elseif (str.eq."4525-111235") then
         include "leshouches_R_200.inc"
         iflow=nint(jamp2200(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2200(i)
         enddo
         goto 20
      elseif (str.eq."4525-111255") then
         include "leshouches_R_201.inc"
         iflow=nint(jamp2201(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2201(i)
         enddo
         goto 20
      elseif (str.eq."4025-111210") then
         include "leshouches_R_202.inc"
         iflow=nint(jamp2202(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2202(i)
         enddo
         goto 20
      elseif (str.eq."4025-111230") then
         include "leshouches_R_203.inc"
         iflow=nint(jamp2203(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2203(i)
         enddo
         goto 20
      elseif (str.eq."4025-111250") then
         include "leshouches_R_204.inc"
         iflow=nint(jamp2204(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2204(i)
         enddo
         goto 20
      elseif (str.eq."-3-125-1112-1-2") then
         include "leshouches_R_205.inc"
         iflow=nint(jamp2205(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2205(i)
         enddo
         goto 20
      elseif (str.eq."-3-125-1112-1-4") then
         include "leshouches_R_206.inc"
         iflow=nint(jamp2206(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2206(i)
         enddo
         goto 20
      elseif (str.eq."-3-125-1112-2-3") then
         include "leshouches_R_207.inc"
         iflow=nint(jamp2207(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2207(i)
         enddo
         goto 20
      elseif (str.eq."-3-125-1112-4-3") then
         include "leshouches_R_208.inc"
         iflow=nint(jamp2208(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2208(i)
         enddo
         goto 20
      elseif (str.eq."-3125-11121-2") then
         include "leshouches_R_209.inc"
         iflow=nint(jamp2209(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2209(i)
         enddo
         goto 20
      elseif (str.eq."-3125-11121-4") then
         include "leshouches_R_210.inc"
         iflow=nint(jamp2210(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2210(i)
         enddo
         goto 20
      elseif (str.eq."-3-225-1112-2-2") then
         include "leshouches_R_211.inc"
         iflow=nint(jamp2211(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2211(i)
         enddo
         goto 20
      elseif (str.eq."-3-225-1112-2-4") then
         include "leshouches_R_212.inc"
         iflow=nint(jamp2212(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2212(i)
         enddo
         goto 20
      elseif (str.eq."-3225-11121-1") then
         include "leshouches_R_213.inc"
         iflow=nint(jamp2213(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2213(i)
         enddo
         goto 20
      elseif (str.eq."-3225-11121-3") then
         include "leshouches_R_214.inc"
         iflow=nint(jamp2214(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2214(i)
         enddo
         goto 20
      elseif (str.eq."-3225-11122-2") then
         include "leshouches_R_215.inc"
         iflow=nint(jamp2215(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2215(i)
         enddo
         goto 20
      elseif (str.eq."-3225-11122-4") then
         include "leshouches_R_216.inc"
         iflow=nint(jamp2216(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2216(i)
         enddo
         goto 20
      elseif (str.eq."-3225-11124-4") then
         include "leshouches_R_217.inc"
         iflow=nint(jamp2217(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2217(i)
         enddo
         goto 20
      elseif (str.eq."-3225-11123-3") then
         include "leshouches_R_218.inc"
         iflow=nint(jamp2218(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2218(i)
         enddo
         goto 20
      elseif (str.eq."-3225-1112-35") then
         include "leshouches_R_219.inc"
         iflow=nint(jamp2219(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2219(i)
         enddo
         goto 20
      elseif (str.eq."-3225-11125-5") then
         include "leshouches_R_220.inc"
         iflow=nint(jamp2220(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2220(i)
         enddo
         goto 20
      elseif (str.eq."-3225-111200") then
         include "leshouches_R_221.inc"
         iflow=nint(jamp2221(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2221(i)
         enddo
         goto 20
      elseif (str.eq."-3-425-1112-2-4") then
         include "leshouches_R_222.inc"
         iflow=nint(jamp2222(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2222(i)
         enddo
         goto 20
      elseif (str.eq."-3-425-1112-4-4") then
         include "leshouches_R_223.inc"
         iflow=nint(jamp2223(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2223(i)
         enddo
         goto 20
      elseif (str.eq."-3425-11121-1") then
         include "leshouches_R_224.inc"
         iflow=nint(jamp2224(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2224(i)
         enddo
         goto 20
      elseif (str.eq."-3425-11121-3") then
         include "leshouches_R_225.inc"
         iflow=nint(jamp2225(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2225(i)
         enddo
         goto 20
      elseif (str.eq."-3425-11122-2") then
         include "leshouches_R_226.inc"
         iflow=nint(jamp2226(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2226(i)
         enddo
         goto 20
      elseif (str.eq."-3425-1112-24") then
         include "leshouches_R_227.inc"
         iflow=nint(jamp2227(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2227(i)
         enddo
         goto 20
      elseif (str.eq."-3425-11124-4") then
         include "leshouches_R_228.inc"
         iflow=nint(jamp2228(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2228(i)
         enddo
         goto 20
      elseif (str.eq."-3425-11123-3") then
         include "leshouches_R_229.inc"
         iflow=nint(jamp2229(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2229(i)
         enddo
         goto 20
      elseif (str.eq."-3425-1112-35") then
         include "leshouches_R_230.inc"
         iflow=nint(jamp2230(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2230(i)
         enddo
         goto 20
      elseif (str.eq."-3425-11125-5") then
         include "leshouches_R_231.inc"
         iflow=nint(jamp2231(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2231(i)
         enddo
         goto 20
      elseif (str.eq."-3425-111200") then
         include "leshouches_R_232.inc"
         iflow=nint(jamp2232(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2232(i)
         enddo
         goto 20
      elseif (str.eq."-3-325-1112-2-3") then
         include "leshouches_R_233.inc"
         iflow=nint(jamp2233(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2233(i)
         enddo
         goto 20
      elseif (str.eq."-3-325-1112-4-3") then
         include "leshouches_R_234.inc"
         iflow=nint(jamp2234(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2234(i)
         enddo
         goto 20
      elseif (str.eq."-3325-11121-2") then
         include "leshouches_R_235.inc"
         iflow=nint(jamp2235(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2235(i)
         enddo
         goto 20
      elseif (str.eq."-3325-11121-4") then
         include "leshouches_R_236.inc"
         iflow=nint(jamp2236(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2236(i)
         enddo
         goto 20
      elseif (str.eq."-3325-1112-23") then
         include "leshouches_R_237.inc"
         iflow=nint(jamp2237(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2237(i)
         enddo
         goto 20
      elseif (str.eq."-3325-1112-25") then
         include "leshouches_R_238.inc"
         iflow=nint(jamp2238(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2238(i)
         enddo
         goto 20
      elseif (str.eq."-3325-1112-43") then
         include "leshouches_R_239.inc"
         iflow=nint(jamp2239(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2239(i)
         enddo
         goto 20
      elseif (str.eq."-3325-1112-45") then
         include "leshouches_R_240.inc"
         iflow=nint(jamp2240(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2240(i)
         enddo
         goto 20
      elseif (str.eq."-3-525-1112-2-3") then
         include "leshouches_R_241.inc"
         iflow=nint(jamp2241(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2241(i)
         enddo
         goto 20
      elseif (str.eq."-3-525-1112-2-5") then
         include "leshouches_R_242.inc"
         iflow=nint(jamp2242(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2242(i)
         enddo
         goto 20
      elseif (str.eq."-3-525-1112-4-3") then
         include "leshouches_R_243.inc"
         iflow=nint(jamp2243(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2243(i)
         enddo
         goto 20
      elseif (str.eq."-3-525-1112-4-5") then
         include "leshouches_R_244.inc"
         iflow=nint(jamp2244(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2244(i)
         enddo
         goto 20
      elseif (str.eq."-3525-1112-25") then
         include "leshouches_R_245.inc"
         iflow=nint(jamp2245(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2245(i)
         enddo
         goto 20
      elseif (str.eq."-3525-1112-45") then
         include "leshouches_R_246.inc"
         iflow=nint(jamp2246(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2246(i)
         enddo
         goto 20
      elseif (str.eq."-3025-1112-20") then
         include "leshouches_R_247.inc"
         iflow=nint(jamp2247(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2247(i)
         enddo
         goto 20
      elseif (str.eq."-3025-1112-40") then
         include "leshouches_R_248.inc"
         iflow=nint(jamp2248(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2248(i)
         enddo
         goto 20
      elseif (str.eq."3-125-1112-23") then
         include "leshouches_R_249.inc"
         iflow=nint(jamp2249(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2249(i)
         enddo
         goto 20
      elseif (str.eq."3-125-1112-43") then
         include "leshouches_R_250.inc"
         iflow=nint(jamp2250(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2250(i)
         enddo
         goto 20
      elseif (str.eq."3225-111213") then
         include "leshouches_R_251.inc"
         iflow=nint(jamp2251(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2251(i)
         enddo
         goto 20
      elseif (str.eq."3225-111233") then
         include "leshouches_R_252.inc"
         iflow=nint(jamp2252(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2252(i)
         enddo
         goto 20
      elseif (str.eq."3225-111235") then
         include "leshouches_R_253.inc"
         iflow=nint(jamp2253(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2253(i)
         enddo
         goto 20
      elseif (str.eq."3425-111213") then
         include "leshouches_R_254.inc"
         iflow=nint(jamp2254(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2254(i)
         enddo
         goto 20
      elseif (str.eq."3425-111233") then
         include "leshouches_R_255.inc"
         iflow=nint(jamp2255(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2255(i)
         enddo
         goto 20
      elseif (str.eq."3425-111235") then
         include "leshouches_R_256.inc"
         iflow=nint(jamp2256(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2256(i)
         enddo
         goto 20
      elseif (str.eq."3-325-11121-2") then
         include "leshouches_R_257.inc"
         iflow=nint(jamp2257(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2257(i)
         enddo
         goto 20
      elseif (str.eq."3-325-11121-4") then
         include "leshouches_R_258.inc"
         iflow=nint(jamp2258(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2258(i)
         enddo
         goto 20
      elseif (str.eq."3-325-1112-23") then
         include "leshouches_R_259.inc"
         iflow=nint(jamp2259(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2259(i)
         enddo
         goto 20
      elseif (str.eq."3-325-1112-25") then
         include "leshouches_R_260.inc"
         iflow=nint(jamp2260(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2260(i)
         enddo
         goto 20
      elseif (str.eq."3-325-1112-43") then
         include "leshouches_R_261.inc"
         iflow=nint(jamp2261(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2261(i)
         enddo
         goto 20
      elseif (str.eq."3-325-1112-45") then
         include "leshouches_R_262.inc"
         iflow=nint(jamp2262(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2262(i)
         enddo
         goto 20
      elseif (str.eq."3-525-1112-23") then
         include "leshouches_R_263.inc"
         iflow=nint(jamp2263(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2263(i)
         enddo
         goto 20
      elseif (str.eq."3-525-1112-43") then
         include "leshouches_R_264.inc"
         iflow=nint(jamp2264(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2264(i)
         enddo
         goto 20
      elseif (str.eq."-5-125-1112-1-2") then
         include "leshouches_R_265.inc"
         iflow=nint(jamp2265(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2265(i)
         enddo
         goto 20
      elseif (str.eq."-5-125-1112-1-4") then
         include "leshouches_R_266.inc"
         iflow=nint(jamp2266(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2266(i)
         enddo
         goto 20
      elseif (str.eq."-5-125-1112-2-5") then
         include "leshouches_R_267.inc"
         iflow=nint(jamp2267(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2267(i)
         enddo
         goto 20
      elseif (str.eq."-5-125-1112-4-5") then
         include "leshouches_R_268.inc"
         iflow=nint(jamp2268(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2268(i)
         enddo
         goto 20
      elseif (str.eq."-5125-11121-2") then
         include "leshouches_R_269.inc"
         iflow=nint(jamp2269(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2269(i)
         enddo
         goto 20
      elseif (str.eq."-5125-11121-4") then
         include "leshouches_R_270.inc"
         iflow=nint(jamp2270(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2270(i)
         enddo
         goto 20
      elseif (str.eq."-5-225-1112-2-2") then
         include "leshouches_R_271.inc"
         iflow=nint(jamp2271(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2271(i)
         enddo
         goto 20
      elseif (str.eq."-5-225-1112-2-4") then
         include "leshouches_R_272.inc"
         iflow=nint(jamp2272(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2272(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11121-1") then
         include "leshouches_R_273.inc"
         iflow=nint(jamp2273(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2273(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11121-5") then
         include "leshouches_R_274.inc"
         iflow=nint(jamp2274(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2274(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11122-2") then
         include "leshouches_R_275.inc"
         iflow=nint(jamp2275(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2275(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11122-4") then
         include "leshouches_R_276.inc"
         iflow=nint(jamp2276(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2276(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11124-4") then
         include "leshouches_R_277.inc"
         iflow=nint(jamp2277(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2277(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11123-3") then
         include "leshouches_R_278.inc"
         iflow=nint(jamp2278(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2278(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11123-5") then
         include "leshouches_R_279.inc"
         iflow=nint(jamp2279(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2279(i)
         enddo
         goto 20
      elseif (str.eq."-5225-11125-5") then
         include "leshouches_R_280.inc"
         iflow=nint(jamp2280(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2280(i)
         enddo
         goto 20
      elseif (str.eq."-5225-111200") then
         include "leshouches_R_281.inc"
         iflow=nint(jamp2281(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2281(i)
         enddo
         goto 20
      elseif (str.eq."-5-425-1112-2-4") then
         include "leshouches_R_282.inc"
         iflow=nint(jamp2282(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2282(i)
         enddo
         goto 20
      elseif (str.eq."-5-425-1112-4-4") then
         include "leshouches_R_283.inc"
         iflow=nint(jamp2283(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2283(i)
         enddo
         goto 20
      elseif (str.eq."-5425-11121-1") then
         include "leshouches_R_284.inc"
         iflow=nint(jamp2284(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2284(i)
         enddo
         goto 20
      elseif (str.eq."-5425-11121-5") then
         include "leshouches_R_285.inc"
         iflow=nint(jamp2285(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2285(i)
         enddo
         goto 20
      elseif (str.eq."-5425-11122-2") then
         include "leshouches_R_286.inc"
         iflow=nint(jamp2286(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2286(i)
         enddo
         goto 20
      elseif (str.eq."-5425-1112-24") then
         include "leshouches_R_287.inc"
         iflow=nint(jamp2287(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2287(i)
         enddo
         goto 20
      elseif (str.eq."-5425-11124-4") then
         include "leshouches_R_288.inc"
         iflow=nint(jamp2288(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2288(i)
         enddo
         goto 20
      elseif (str.eq."-5425-11123-3") then
         include "leshouches_R_289.inc"
         iflow=nint(jamp2289(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2289(i)
         enddo
         goto 20
      elseif (str.eq."-5425-11123-5") then
         include "leshouches_R_290.inc"
         iflow=nint(jamp2290(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2290(i)
         enddo
         goto 20
      elseif (str.eq."-5425-11125-5") then
         include "leshouches_R_291.inc"
         iflow=nint(jamp2291(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2291(i)
         enddo
         goto 20
      elseif (str.eq."-5425-111200") then
         include "leshouches_R_292.inc"
         iflow=nint(jamp2292(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2292(i)
         enddo
         goto 20
      elseif (str.eq."-5-325-1112-2-3") then
         include "leshouches_R_293.inc"
         iflow=nint(jamp2293(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2293(i)
         enddo
         goto 20
      elseif (str.eq."-5-325-1112-2-5") then
         include "leshouches_R_294.inc"
         iflow=nint(jamp2294(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2294(i)
         enddo
         goto 20
      elseif (str.eq."-5-325-1112-4-3") then
         include "leshouches_R_295.inc"
         iflow=nint(jamp2295(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2295(i)
         enddo
         goto 20
      elseif (str.eq."-5-325-1112-4-5") then
         include "leshouches_R_296.inc"
         iflow=nint(jamp2296(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2296(i)
         enddo
         goto 20
      elseif (str.eq."-5325-1112-23") then
         include "leshouches_R_297.inc"
         iflow=nint(jamp2297(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2297(i)
         enddo
         goto 20
      elseif (str.eq."-5325-1112-43") then
         include "leshouches_R_298.inc"
         iflow=nint(jamp2298(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2298(i)
         enddo
         goto 20
      elseif (str.eq."-5-525-1112-2-5") then
         include "leshouches_R_299.inc"
         iflow=nint(jamp2299(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2299(i)
         enddo
         goto 20
      elseif (str.eq."-5-525-1112-4-5") then
         include "leshouches_R_300.inc"
         iflow=nint(jamp2300(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2300(i)
         enddo
         goto 20
      elseif (str.eq."-5525-11121-2") then
         include "leshouches_R_301.inc"
         iflow=nint(jamp2301(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2301(i)
         enddo
         goto 20
      elseif (str.eq."-5525-11121-4") then
         include "leshouches_R_302.inc"
         iflow=nint(jamp2302(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2302(i)
         enddo
         goto 20
      elseif (str.eq."-5525-1112-23") then
         include "leshouches_R_303.inc"
         iflow=nint(jamp2303(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2303(i)
         enddo
         goto 20
      elseif (str.eq."-5525-1112-25") then
         include "leshouches_R_304.inc"
         iflow=nint(jamp2304(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2304(i)
         enddo
         goto 20
      elseif (str.eq."-5525-1112-43") then
         include "leshouches_R_305.inc"
         iflow=nint(jamp2305(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2305(i)
         enddo
         goto 20
      elseif (str.eq."-5525-1112-45") then
         include "leshouches_R_306.inc"
         iflow=nint(jamp2306(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2306(i)
         enddo
         goto 20
      elseif (str.eq."-5025-1112-20") then
         include "leshouches_R_307.inc"
         iflow=nint(jamp2307(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2307(i)
         enddo
         goto 20
      elseif (str.eq."-5025-1112-40") then
         include "leshouches_R_308.inc"
         iflow=nint(jamp2308(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2308(i)
         enddo
         goto 20
      elseif (str.eq."5-125-1112-25") then
         include "leshouches_R_309.inc"
         iflow=nint(jamp2309(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2309(i)
         enddo
         goto 20
      elseif (str.eq."5-125-1112-45") then
         include "leshouches_R_310.inc"
         iflow=nint(jamp2310(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2310(i)
         enddo
         goto 20
      elseif (str.eq."5225-111215") then
         include "leshouches_R_311.inc"
         iflow=nint(jamp2311(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2311(i)
         enddo
         goto 20
      elseif (str.eq."5225-111235") then
         include "leshouches_R_312.inc"
         iflow=nint(jamp2312(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2312(i)
         enddo
         goto 20
      elseif (str.eq."5225-111255") then
         include "leshouches_R_313.inc"
         iflow=nint(jamp2313(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2313(i)
         enddo
         goto 20
      elseif (str.eq."5425-111215") then
         include "leshouches_R_314.inc"
         iflow=nint(jamp2314(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2314(i)
         enddo
         goto 20
      elseif (str.eq."5425-111235") then
         include "leshouches_R_315.inc"
         iflow=nint(jamp2315(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2315(i)
         enddo
         goto 20
      elseif (str.eq."5425-111255") then
         include "leshouches_R_316.inc"
         iflow=nint(jamp2316(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2316(i)
         enddo
         goto 20
      elseif (str.eq."5-325-1112-25") then
         include "leshouches_R_317.inc"
         iflow=nint(jamp2317(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2317(i)
         enddo
         goto 20
      elseif (str.eq."5-325-1112-45") then
         include "leshouches_R_318.inc"
         iflow=nint(jamp2318(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2318(i)
         enddo
         goto 20
      elseif (str.eq."5-525-11121-2") then
         include "leshouches_R_319.inc"
         iflow=nint(jamp2319(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2319(i)
         enddo
         goto 20
      elseif (str.eq."5-525-11121-4") then
         include "leshouches_R_320.inc"
         iflow=nint(jamp2320(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2320(i)
         enddo
         goto 20
      elseif (str.eq."5-525-1112-23") then
         include "leshouches_R_321.inc"
         iflow=nint(jamp2321(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2321(i)
         enddo
         goto 20
      elseif (str.eq."5-525-1112-25") then
         include "leshouches_R_322.inc"
         iflow=nint(jamp2322(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2322(i)
         enddo
         goto 20
      elseif (str.eq."5-525-1112-43") then
         include "leshouches_R_323.inc"
         iflow=nint(jamp2323(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2323(i)
         enddo
         goto 20
      elseif (str.eq."5-525-1112-45") then
         include "leshouches_R_324.inc"
         iflow=nint(jamp2324(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2324(i)
         enddo
         goto 20
      elseif (str.eq."0-125-1112-20") then
         include "leshouches_R_325.inc"
         iflow=nint(jamp2325(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2325(i)
         enddo
         goto 20
      elseif (str.eq."0-125-1112-40") then
         include "leshouches_R_326.inc"
         iflow=nint(jamp2326(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2326(i)
         enddo
         goto 20
      elseif (str.eq."0225-111210") then
         include "leshouches_R_327.inc"
         iflow=nint(jamp2327(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2327(i)
         enddo
         goto 20
      elseif (str.eq."0225-111230") then
         include "leshouches_R_328.inc"
         iflow=nint(jamp2328(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2328(i)
         enddo
         goto 20
      elseif (str.eq."0225-111250") then
         include "leshouches_R_329.inc"
         iflow=nint(jamp2329(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2329(i)
         enddo
         goto 20
      elseif (str.eq."0425-111210") then
         include "leshouches_R_330.inc"
         iflow=nint(jamp2330(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2330(i)
         enddo
         goto 20
      elseif (str.eq."0425-111230") then
         include "leshouches_R_331.inc"
         iflow=nint(jamp2331(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2331(i)
         enddo
         goto 20
      elseif (str.eq."0425-111250") then
         include "leshouches_R_332.inc"
         iflow=nint(jamp2332(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2332(i)
         enddo
         goto 20
      elseif (str.eq."0-325-1112-20") then
         include "leshouches_R_333.inc"
         iflow=nint(jamp2333(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2333(i)
         enddo
         goto 20
      elseif (str.eq."0-325-1112-40") then
         include "leshouches_R_334.inc"
         iflow=nint(jamp2334(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2334(i)
         enddo
         goto 20
      elseif (str.eq."0-525-1112-20") then
         include "leshouches_R_335.inc"
         iflow=nint(jamp2335(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2335(i)
         enddo
         goto 20
      elseif (str.eq."0-525-1112-40") then
         include "leshouches_R_336.inc"
         iflow=nint(jamp2336(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2336(i)
         enddo
         goto 20
      elseif (str.eq."0025-11121-2") then
         include "leshouches_R_337.inc"
         iflow=nint(jamp2337(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2337(i)
         enddo
         goto 20
      elseif (str.eq."0025-11121-4") then
         include "leshouches_R_338.inc"
         iflow=nint(jamp2338(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2338(i)
         enddo
         goto 20
      elseif (str.eq."0025-1112-23") then
         include "leshouches_R_339.inc"
         iflow=nint(jamp2339(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2339(i)
         enddo
         goto 20
      elseif (str.eq."0025-1112-25") then
         include "leshouches_R_340.inc"
         iflow=nint(jamp2340(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2340(i)
         enddo
         goto 20
      elseif (str.eq."0025-1112-43") then
         include "leshouches_R_341.inc"
         iflow=nint(jamp2341(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2341(i)
         enddo
         goto 20
      elseif (str.eq."0025-1112-45") then
         include "leshouches_R_342.inc"
         iflow=nint(jamp2342(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2342(i)
         enddo
         goto 20
      endif
      
      do while(mtc)
         do i= 5+1,nexternal
            ic(i)=ic(i)- 5
         enddo
         goto 10
      enddo
      if (.not.mtc) then
         write (*,*) "Error #1, in sborn_proc.f"
         stop
      endif
      
 20   continue
      xtarget=jamp2cum(iflow)*random()
      ifl=1
      do while (jamp2cum(ifl).lt.xtarget)
         ifl=ifl+1
      enddo
      do i=1,2
         do j=1,nexternal
            color1(i,j)=ICOLUP(i,j,ifl)
         enddo
      enddo
      call switchcolor(color1,color,
     &     ic,nexternal)
      
      return
      end
      
      
      
      
