      f[1][1] = 0.0;
      f[1][2] = 0.0;
      f[1][3] = 0.0;
      f[1][4] = 0.0;
      f[1][5] = 0.0;
      f[1][6] = 0.0;
      f[1][7] = 0.0;
      f[1][8] = 0.0;
      f[1][9] = 0.0;
      t1 = Bptr[3][B1];
      t2 = m*m;
      t4 = dotp(kq,kqb);
      t6 = Bptr[3][B00];
      t8 = nfl*(t1*t2+t1*t4-t6);
      t9 = dotp(p2,p3);
      t10 = 1.0/t9;
      t11 = t2+t4;
      t12 = t11*t11;
      t13 = 1.0/t12;
      t15 = t8*t10*t13;
      t16 = t15/4.0;
      f[1][10] = -t16;
      f[1][11] = 0.0;
      f[1][12] = 0.0;
      f[1][13] = 0.0;
      f[1][14] = 0.0;
      f[1][15] = 0.0;
      f[1][16] = 0.0;
      f[1][17] = 0.0;
      f[1][18] = 0.0;
      f[1][19] = 0.0;
      f[1][20] = 0.0;
      f[1][21] = 0.0;
      f[1][22] = 0.0;
      f[1][23] = 0.0;
      f[1][24] = 0.0;
      t17 = dotp(p1,p3);
      t20 = 1.0/t17;
      t21 = t20*t13;
      f[1][25] = t8*(t9+t17)*t21*t10/8.0;
      f[1][26] = 0.0;
      f[1][27] = 0.0;
      f[1][28] = 0.0;
      f[1][29] = 0.0;
      f[1][30] = 0.0;
      f[1][31] = 0.0;
      f[1][32] = 0.0;
      f[1][33] = 0.0;
      f[1][34] = 0.0;
      f[1][35] = 0.0;
      f[1][36] = 0.0;
      f[1][37] = 0.0;
      f[1][38] = 0.0;
      f[1][39] = 0.0;
      f[1][40] = 0.0;
      f[1][41] = 0.0;
      f[1][42] = 0.0;
      f[1][43] = 0.0;
      f[1][44] = 0.0;
      f[1][45] = 0.0;
      f[1][46] = f[1][10];
      f[1][47] = 0.0;
      f[1][48] = 0.0;
      f[1][49] = 0.0;
      f[1][50] = f[1][25];
      f[1][51] = 0.0;
      f[1][52] = 0.0;
      f[1][53] = 0.0;
      f[1][54] = 0.0;
      f[1][55] = 0.0;
      f[1][56] = 0.0;
      f[1][57] = 0.0;
      f[1][58] = 0.0;
      f[1][59] = 0.0;
      f[1][60] = 0.0;
      f[1][61] = 0.0;
      f[1][62] = 0.0;
      f[1][63] = f[1][46];
      t24 = t8*t21;
      f[1][64] = t24/4.0;
      f[1][65] = 0.0;
      f[1][66] = 0.0;
      f[1][67] = 0.0;
      f[1][68] = 0.0;
      f[1][69] = 0.0;
      f[1][70] = 0.0;
      f[1][71] = 0.0;
      f[1][72] = 0.0;
      f[1][73] = 0.0;
      f[1][74] = 0.0;
      f[1][75] = t16;
      f[1][76] = f[1][75];
      f[1][77] = 0.0;
      f[1][78] = 0.0;
      f[1][79] = 0.0;
      f[1][80] = 0.0;
      f[1][81] = 0.0;
      f[1][82] = 0.0;
      f[1][83] = 0.0;
      f[1][84] = 0.0;
      f[1][85] = 0.0;
      f[2][1] = 0.0;
      t25 = Cptr[8][C22];
      t26 = Cptr[8][C12];
      t27 = Cptr[8][C112];
      t28 = Cptr[8][C122];
      f[2][2] = nfl*(t25-t26-t27+t28)/t11/s;
      f[2][3] = f[2][2];
      f[2][4] = 0.0;
      f[2][5] = f[2][3];
      f[2][6] = f[2][5];
      f[2][7] = 0.0;
      f[2][8] = 0.0;
      f[2][9] = 0.0;
      t34 = Bptr[1][B00];
      t35 = t34*t9;
      t36 = t2*t2;
      t39 = s*s;
      t40 = t6*t39;
      t41 = dotp(p3,kqb);
      t43 = t17*t26;
      t44 = t43*s;
      t45 = t9*t41;
      t46 = t45*t2;
      t49 = t4*t4;
      t52 = t41*t2;
      t55 = t41*t4;
      t58 = t1*t39;
      t61 = Bptr[1][B1];
      t62 = t61*s;
      t69 = t2*t4;
      t72 = t6*s;
      t75 = t17*t25;
      t76 = t75*s;
      t79 = t45*t4;
      t84 = t9*t9;
      t85 = t84*t25;
      t86 = s*t41;
      t87 = t86*t2;
      t90 = t86*t4;
      t93 = -4.0*t35*t36+t40*t41-2.0*t44*t46-4.0*t35*t49-4.0*t35*t52-4.0*t35*
t55-t58*t52-t58*t55+2.0*t62*t9*t36+2.0*t62*t9*t49-8.0*t35*t69-2.0*t72*t45+2.0*
t76*t46+2.0*t76*t79-2.0*t44*t79+2.0*t85*t87+2.0*t85*t90;
      t94 = t84*t26;
      t103 = Bptr[3][B0];
      t104 = t103*s;
      t109 = t1*s;
      t114 = Cptr[8][C2];
      t115 = t114*t39;
      t118 = Cptr[8][C00];
      t119 = t118*s;
      t124 = Cptr[8][C001];
      t125 = t124*s;
      t130 = Cptr[8][C002];
      t131 = t130*s;
      t136 = t9*t2;
      t140 = -2.0*t94*t87-2.0*t94*t90+2.0*t62*t46+2.0*t62*t79+2.0*t104*t46+2.0*
t104*t79+3.0*t109*t46+3.0*t109*t79+t115*t46+t115*t79-2.0*t119*t46-2.0*t119*t79+
4.0*t125*t46+4.0*t125*t79-4.0*t131*t46-4.0*t131*t79+4.0*t62*t136*t4;
      t144 = 1.0/t39;
      t145 = t144*t10;
      t146 = 1.0/t41;
      f[2][10] = -nfl*(t93+t140)*t13*t145*t146/4.0;
      f[2][11] = 0.0;
      f[2][12] = 0.0;
      f[2][13] = 0.0;
      f[2][14] = 0.0;
      f[2][15] = 0.0;
      f[2][16] = 0.0;
      f[2][17] = 0.0;
      f[2][18] = 0.0;
      f[2][19] = 0.0;
      f[2][20] = 0.0;
      f[2][21] = 0.0;
      f[2][22] = 0.0;
      f[2][23] = 0.0;
      f[2][24] = 0.0;
      f[2][25] = -t15/8.0;
      f[2][26] = 0.0;
      f[2][27] = 0.0;
      f[2][28] = 0.0;
      f[2][29] = 0.0;
      f[2][30] = 0.0;
      f[2][31] = 0.0;
      f[2][32] = 0.0;
      f[2][33] = 0.0;
      f[2][34] = 0.0;
      f[2][35] = 0.0;
      f[2][36] = 0.0;
      f[2][37] = 0.0;
      f[2][38] = f[2][6];
      f[2][39] = f[2][38];
      f[2][40] = 0.0;
      f[2][41] = f[2][39];
      f[2][42] = f[2][41];
      f[2][43] = 0.0;
      f[2][44] = 0.0;
      f[2][45] = 0.0;
      t152 = 4.0*t35*t4;
      t154 = 4.0*t35*t2;
      t155 = s*t9;
      t156 = t155*t2;
      t158 = 2.0*t43*t156;
      t159 = t58*t4;
      t161 = 2.0*t72*t9;
      t162 = t58*t2;
      t163 = s*t2;
      t165 = 2.0*t85*t163;
      t167 = 2.0*t62*t136;
      t168 = t9*t4;
      t170 = 2.0*t62*t168;
      t171 = s*t4;
      t173 = 2.0*t85*t171;
      t175 = 2.0*t94*t163;
      t177 = 2.0*t94*t171;
      t178 = t115*t136;
      t179 = -t152-t154-t158+t40-t159-t161-t162+t165+t167+t170+t173-t175-t177+
t178;
      t180 = t115*t168;
      t181 = t119*t136;
      t182 = 2.0*t181;
      t183 = t104*t136;
      t185 = t104*t168;
      t187 = t109*t136;
      t189 = t109*t168;
      t191 = t119*t168;
      t192 = 2.0*t191;
      t194 = 4.0*t125*t136;
      t196 = 4.0*t125*t168;
      t198 = 4.0*t131*t136;
      t200 = 4.0*t131*t168;
      t202 = 2.0*t75*t156;
      t203 = t155*t4;
      t205 = 2.0*t75*t203;
      t207 = 2.0*t43*t203;
      t208 = t180-t182+2.0*t183+2.0*t185+3.0*t187+3.0*t189-t192+t194+t196-t198-
t200+t202+t205-t207;
      t211 = t145*t13;
      f[2][46] = -nfl*(t179+t208)*t211/4.0;
      f[2][47] = 0.0;
      f[2][48] = 0.0;
      f[2][49] = 0.0;
      f[2][50] = f[2][25];
      f[2][51] = 0.0;
      f[2][52] = 0.0;
      f[2][53] = 0.0;
      f[2][54] = 0.0;
      f[2][55] = 0.0;
      f[2][56] = 0.0;
      f[2][57] = 0.0;
      f[2][58] = 0.0;
      f[2][59] = 0.0;
      f[2][60] = 0.0;
      f[2][61] = 0.0;
      f[2][62] = 0.0;
      t218 = t167+t170-t154-t152+t183+t185+2.0*t187+2.0*t189+t178+t180-4.0*t181
-4.0*t191-t194-t196-t161-t162-t159+t40;
      f[2][63] = -nfl*t218*t211/4.0;
      t223 = 2.0*t62*t2;
      t225 = 2.0*t62*t4;
      t226 = t34*t2;
      t227 = 4.0*t226;
      t229 = 4.0*t34*t4;
      t230 = t104*t2;
      t231 = t104*t4;
      t232 = t109*t2;
      t234 = t109*t4;
      t236 = t115*t2;
      t237 = t115*t4;
      t238 = t119*t2;
      t240 = t119*t4;
      t243 = 4.0*t125*t2;
      t245 = 4.0*t125*t4;
      t246 = 2.0*t72;
      t247 = t223+t225-t227-t229+t230+t231+2.0*t232+2.0*t234+t236+t237-4.0*t238
-4.0*t240-t243-t245-t246;
      t249 = t144*t13;
      f[2][64] = -nfl*t247*t249/4.0;
      f[2][65] = 0.0;
      f[2][66] = 0.0;
      f[2][67] = 0.0;
      f[2][68] = 0.0;
      f[2][69] = 0.0;
      f[2][70] = 0.0;
      t253 = t62-2.0*t34;
      t254 = nfl*t253;
      t256 = t254*t146*t144;
      f[2][71] = -t256/4.0;
      f[2][72] = 0.0;
      f[2][73] = 0.0;
      f[2][74] = f[2][71];
      t258 = t152+t154+t158-t40+t159+t161+t162-t165-t167-t170-t173+t175+t177;
      t259 = t25*t9;
      t266 = t182-t183-t185-t187-t189+t192-t198-t200-t202-t205+t207+2.0*t259*
t39*t2+2.0*t259*t39*t4;
      f[2][75] = -nfl*(t258+t266)*t211/4.0;
      f[2][76] = f[2][75];
      f[2][77] = 0.0;
      f[2][78] = 0.0;
      f[2][79] = 0.0;
      f[2][80] = 0.0;
      f[2][81] = 0.0;
      f[2][82] = 0.0;
      f[2][83] = 0.0;
      f[2][84] = 0.0;
      f[2][85] = 0.0;
      f[3][1] = 0.0;
      f[3][2] = 0.0;
      f[3][3] = 0.0;
      f[3][4] = 0.0;
      f[3][5] = 0.0;
      f[3][6] = 0.0;
      f[3][7] = 0.0;
      f[3][8] = 0.0;
      f[3][9] = 0.0;
      f[3][10] = t256/2.0;
      f[3][11] = 0.0;
      f[3][12] = 0.0;
      f[3][13] = 0.0;
      f[3][14] = 0.0;
      f[3][15] = 0.0;
      f[3][16] = 0.0;
      f[3][17] = 0.0;
      f[3][18] = 0.0;
      f[3][19] = 0.0;
      f[3][20] = 0.0;
      f[3][21] = 0.0;
      f[3][22] = 0.0;
      f[3][23] = 0.0;
      f[3][24] = 0.0;
      f[3][25] = 0.0;
      f[3][26] = 0.0;
      f[3][27] = 0.0;
      f[3][28] = 0.0;
      f[3][29] = 0.0;
      f[3][30] = 0.0;
      f[3][31] = 0.0;
      f[3][32] = 0.0;
      f[3][33] = 0.0;
      f[3][34] = 0.0;
      f[3][35] = 0.0;
      f[3][36] = 0.0;
      f[3][37] = 0.0;
      f[3][38] = 0.0;
      f[3][39] = 0.0;
      f[3][40] = 0.0;
      f[3][41] = 0.0;
      f[3][42] = 0.0;
      f[3][43] = 0.0;
      f[3][44] = 0.0;
      f[3][45] = 0.0;
      t271 = dotp(p3,kq);
      t272 = 1.0/t271;
      t273 = t272*t144;
      t274 = t254*t273;
      t275 = t274/2.0;
      f[3][46] = -t275;
      f[3][47] = 0.0;
      f[3][48] = 0.0;
      f[3][49] = 0.0;
      f[3][50] = 0.0;
      f[3][51] = 0.0;
      f[3][52] = 0.0;
      f[3][53] = 0.0;
      f[3][54] = 0.0;
      f[3][55] = 0.0;
      f[3][56] = 0.0;
      f[3][57] = 0.0;
      f[3][58] = 0.0;
      f[3][59] = 0.0;
      f[3][60] = 0.0;
      f[3][61] = 0.0;
      f[3][62] = 0.0;
      f[3][63] = f[3][46];
      f[3][64] = f[3][63];
      f[3][65] = 0.0;
      f[3][66] = 0.0;
      f[3][67] = 0.0;
      f[3][68] = 0.0;
      f[3][69] = 0.0;
      f[3][70] = 0.0;
      f[3][71] = nfl*(t41+t271)*t253*t273*t146/4.0;
      f[3][72] = 0.0;
      f[3][73] = 0.0;
      f[3][74] = f[3][71];
      f[3][75] = t275;
      f[3][76] = f[3][75];
      f[3][77] = 0.0;
      f[3][78] = 0.0;
      f[3][79] = 0.0;
      f[3][80] = 0.0;
      f[3][81] = 0.0;
      f[3][82] = 0.0;
      f[3][83] = 0.0;
      f[3][84] = 0.0;
      f[3][85] = 0.0;
      f[4][1] = 0.0;
      f[4][2] = -f[2][42];
      f[4][3] = f[4][2];
      f[4][4] = 0.0;
      f[4][5] = f[4][3];
      f[4][6] = f[4][5];
      f[4][7] = 0.0;
      f[4][8] = 0.0;
      f[4][9] = 0.0;
      t287 = t26*s;
      t290 = t25*s;
      t295 = t223-4.0*t131*t2+t236+t237+3.0*t232+t245-2.0*t238-4.0*t131*t4+t243
-2.0*t287*t136+2.0*t290*t168+2.0*t290*t136;
      t310 = -2.0*t287*t168+2.0*t75*t163+2.0*t231+2.0*t75*t171+2.0*t230-t227
-2.0*t43*t163-2.0*t43*t171+t225+3.0*t234-2.0*t240-t246-t229;
      f[4][10] = nfl*(t295+t310)*t249/4.0;
      f[4][11] = 0.0;
      f[4][12] = 0.0;
      f[4][13] = 0.0;
      f[4][14] = 0.0;
      f[4][15] = 0.0;
      f[4][16] = 0.0;
      f[4][17] = 0.0;
      f[4][18] = 0.0;
      f[4][19] = 0.0;
      f[4][20] = 0.0;
      f[4][21] = 0.0;
      f[4][22] = 0.0;
      f[4][23] = 0.0;
      f[4][24] = 0.0;
      f[4][25] = -t24/8.0;
      f[4][26] = 0.0;
      f[4][27] = 0.0;
      f[4][28] = 0.0;
      f[4][29] = 0.0;
      f[4][30] = 0.0;
      f[4][31] = 0.0;
      f[4][32] = 0.0;
      f[4][33] = 0.0;
      f[4][34] = 0.0;
      f[4][35] = 0.0;
      f[4][36] = 0.0;
      f[4][37] = 0.0;
      f[4][38] = f[4][6];
      f[4][39] = f[4][38];
      f[4][40] = 0.0;
      f[4][41] = f[4][39];
      f[4][42] = f[4][41];
      f[4][43] = 0.0;
      f[4][44] = 0.0;
      f[4][45] = 0.0;
      t316 = 2.0*t72*t271;
      t318 = 8.0*t226*t4;
      t319 = t61*t271;
      t321 = 2.0*t319*t171;
      t322 = t34*t271;
      t324 = 4.0*t322*t4;
      t326 = 2.0*t62*t49;
      t327 = t271*t2;
      t328 = t104*t327;
      t330 = t271*t4;
      t331 = t104*t330;
      t333 = t109*t327;
      t335 = t109*t330;
      t337 = t115*t327;
      t338 = t115*t330;
      t340 = 4.0*t322*t2;
      t342 = 4.0*t34*t36;
      t344 = 4.0*t34*t49;
      t345 = t119*t327;
      t346 = 2.0*t345;
      t347 = -t316-t318+t321-t324+t326+2.0*t328+2.0*t331+3.0*t333+3.0*t335+t337
+t338-t340-t342-t344-t346;
      t348 = t119*t330;
      t349 = 2.0*t348;
      t351 = 4.0*t125*t327;
      t353 = 4.0*t125*t330;
      t355 = 4.0*t131*t327;
      t357 = 4.0*t131*t330;
      t359 = 4.0*t62*t69;
      t360 = s*t271;
      t361 = t360*t2;
      t363 = 2.0*t75*t361;
      t365 = 2.0*t319*t163;
      t366 = t9*t26;
      t367 = t360*t4;
      t369 = 2.0*t366*t367;
      t371 = 2.0*t366*t361;
      t373 = 2.0*t259*t367;
      t375 = 2.0*t259*t361;
      t377 = 2.0*t43*t367;
      t379 = 2.0*t43*t361;
      t381 = 2.0*t75*t367;
      t383 = 2.0*t62*t36;
      t384 = -t349+t351+t353-t355-t357+t359+t363+t365-t369-t371+t373+t375-t377-
t379+t381+t383;
      t387 = t249*t272;
      f[4][46] = nfl*(t347+t384)*t387/4.0;
      f[4][47] = 0.0;
      f[4][48] = 0.0;
      f[4][49] = 0.0;
      f[4][50] = f[4][25];
      f[4][51] = 0.0;
      f[4][52] = 0.0;
      f[4][53] = 0.0;
      f[4][54] = 0.0;
      f[4][55] = 0.0;
      f[4][56] = 0.0;
      f[4][57] = 0.0;
      f[4][58] = 0.0;
      f[4][59] = 0.0;
      f[4][60] = 0.0;
      f[4][61] = 0.0;
      f[4][62] = 0.0;
      t394 = -4.0*t345-4.0*t348-t351-t353-t316+t383+t359+t326-t342-t318-t344;
      f[4][63] = nfl*(t365+t321-t340-t324+t328+t331+2.0*t333+2.0*t335+t337+t338+
t394)*t387/4.0;
      t398 = t17*t271;
      t399 = t398*t2;
      t402 = t398*t4;
      t405 = t34*t17;
      t422 = -2.0*t62*t399-2.0*t62*t402+4.0*t405*t327+4.0*t405*t330-t104*t399-
t104*t402-2.0*t109*t399-2.0*t109*t402-t115*t399-t115*t402+4.0*t119*t399+4.0*
t119*t402;
      t448 = 4.0*t125*t399+4.0*t125*t402+2.0*t72*t398+t58*t327+t58*t330-t40*
t271-2.0*t62*t17*t36-4.0*t62*t17*t2*t4-2.0*t62*t17*t49+4.0*t405*t36+8.0*t405*
t69+4.0*t405*t49;
      f[4][64] = -nfl*(t422+t448)*t13*t144*t20*t272/4.0;
      f[4][65] = 0.0;
      f[4][66] = 0.0;
      f[4][67] = 0.0;
      f[4][68] = 0.0;
      f[4][69] = 0.0;
      f[4][70] = 0.0;
      f[4][71] = -t274/4.0;
      f[4][72] = 0.0;
      f[4][73] = 0.0;
      f[4][74] = f[4][71];
      t457 = t316+t318-t321+t324-t326-t328-t331-t333-t335+t340+t342+t344+t346+
t349;
      t458 = t25*t39;
      t463 = -t355-t357-t359-t363-t365+t369+t371-t373-t375+t377+t379-t381-t383+
2.0*t458*t327+2.0*t458*t330;
      f[4][75] = nfl*(t457+t463)*t387/4.0;
      f[4][76] = f[4][75];
      f[4][77] = 0.0;
      f[4][78] = 0.0;
      f[4][79] = 0.0;
      f[4][80] = 0.0;
      f[4][81] = 0.0;
      f[4][82] = 0.0;
      f[4][83] = 0.0;
      f[4][84] = 0.0;
      f[4][85] = 0.0;
