/*
 *    This file was auto-generated using the ACADO Toolkit.
 *    
 *    While ACADO Toolkit is free software released under the terms of
 *    the GNU Lesser General Public License (LGPL), the generated code
 *    as such remains the property of the user who used ACADO Toolkit
 *    to generate this code. In particular, user dependent data of the code
 *    do not inherit the GNU LGPL license. On the other hand, parts of the
 *    generated code that are a direct copy of source code from the
 *    ACADO Toolkit or the software tools it is based on, remain, as derived
 *    work, automatically covered by the LGPL license.
 *    
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *    
 */


#include "nmpc_common.h"




/******************************************************************************/
/*                                                                            */
/* ACADO code generation                                                      */
/*                                                                            */
/******************************************************************************/


int nmpc_modelSimulation(  )
{
int ret;

int lRun1;
int lRun2;
ret = 0;
for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
nmpcWorkspace.state[0] = nmpcVariables.x[lRun1 * 14];
nmpcWorkspace.state[1] = nmpcVariables.x[lRun1 * 14 + 1];
nmpcWorkspace.state[2] = nmpcVariables.x[lRun1 * 14 + 2];
nmpcWorkspace.state[3] = nmpcVariables.x[lRun1 * 14 + 3];
nmpcWorkspace.state[4] = nmpcVariables.x[lRun1 * 14 + 4];
nmpcWorkspace.state[5] = nmpcVariables.x[lRun1 * 14 + 5];
nmpcWorkspace.state[6] = nmpcVariables.x[lRun1 * 14 + 6];
nmpcWorkspace.state[7] = nmpcVariables.x[lRun1 * 14 + 7];
nmpcWorkspace.state[8] = nmpcVariables.x[lRun1 * 14 + 8];
nmpcWorkspace.state[9] = nmpcVariables.x[lRun1 * 14 + 9];
nmpcWorkspace.state[10] = nmpcVariables.x[lRun1 * 14 + 10];
nmpcWorkspace.state[11] = nmpcVariables.x[lRun1 * 14 + 11];
nmpcWorkspace.state[12] = nmpcVariables.x[lRun1 * 14 + 12];
nmpcWorkspace.state[13] = nmpcVariables.x[lRun1 * 14 + 13];

nmpcWorkspace.state[266] = nmpcVariables.u[lRun1 * 4];
nmpcWorkspace.state[267] = nmpcVariables.u[lRun1 * 4 + 1];
nmpcWorkspace.state[268] = nmpcVariables.u[lRun1 * 4 + 2];
nmpcWorkspace.state[269] = nmpcVariables.u[lRun1 * 4 + 3];
nmpcWorkspace.state[270] = nmpcVariables.od[lRun1 * 6];
nmpcWorkspace.state[271] = nmpcVariables.od[lRun1 * 6 + 1];
nmpcWorkspace.state[272] = nmpcVariables.od[lRun1 * 6 + 2];
nmpcWorkspace.state[273] = nmpcVariables.od[lRun1 * 6 + 3];
nmpcWorkspace.state[274] = nmpcVariables.od[lRun1 * 6 + 4];
nmpcWorkspace.state[275] = nmpcVariables.od[lRun1 * 6 + 5];

ret = nmpc_integrate(nmpcWorkspace.state, 1);

nmpcWorkspace.d[lRun1 * 14] = nmpcWorkspace.state[0] - nmpcVariables.x[lRun1 * 14 + 14];
nmpcWorkspace.d[lRun1 * 14 + 1] = nmpcWorkspace.state[1] - nmpcVariables.x[lRun1 * 14 + 15];
nmpcWorkspace.d[lRun1 * 14 + 2] = nmpcWorkspace.state[2] - nmpcVariables.x[lRun1 * 14 + 16];
nmpcWorkspace.d[lRun1 * 14 + 3] = nmpcWorkspace.state[3] - nmpcVariables.x[lRun1 * 14 + 17];
nmpcWorkspace.d[lRun1 * 14 + 4] = nmpcWorkspace.state[4] - nmpcVariables.x[lRun1 * 14 + 18];
nmpcWorkspace.d[lRun1 * 14 + 5] = nmpcWorkspace.state[5] - nmpcVariables.x[lRun1 * 14 + 19];
nmpcWorkspace.d[lRun1 * 14 + 6] = nmpcWorkspace.state[6] - nmpcVariables.x[lRun1 * 14 + 20];
nmpcWorkspace.d[lRun1 * 14 + 7] = nmpcWorkspace.state[7] - nmpcVariables.x[lRun1 * 14 + 21];
nmpcWorkspace.d[lRun1 * 14 + 8] = nmpcWorkspace.state[8] - nmpcVariables.x[lRun1 * 14 + 22];
nmpcWorkspace.d[lRun1 * 14 + 9] = nmpcWorkspace.state[9] - nmpcVariables.x[lRun1 * 14 + 23];
nmpcWorkspace.d[lRun1 * 14 + 10] = nmpcWorkspace.state[10] - nmpcVariables.x[lRun1 * 14 + 24];
nmpcWorkspace.d[lRun1 * 14 + 11] = nmpcWorkspace.state[11] - nmpcVariables.x[lRun1 * 14 + 25];
nmpcWorkspace.d[lRun1 * 14 + 12] = nmpcWorkspace.state[12] - nmpcVariables.x[lRun1 * 14 + 26];
nmpcWorkspace.d[lRun1 * 14 + 13] = nmpcWorkspace.state[13] - nmpcVariables.x[lRun1 * 14 + 27];

for (lRun2 = 0; lRun2 < 196; ++lRun2)
nmpcWorkspace.evGx[(0) + ((lRun2) + (lRun1 * 196))] = nmpcWorkspace.state[lRun2 + 14];


nmpcWorkspace.evGu[lRun1 * 56] = nmpcWorkspace.state[210];
nmpcWorkspace.evGu[lRun1 * 56 + 1] = nmpcWorkspace.state[211];
nmpcWorkspace.evGu[lRun1 * 56 + 2] = nmpcWorkspace.state[212];
nmpcWorkspace.evGu[lRun1 * 56 + 3] = nmpcWorkspace.state[213];
nmpcWorkspace.evGu[lRun1 * 56 + 4] = nmpcWorkspace.state[214];
nmpcWorkspace.evGu[lRun1 * 56 + 5] = nmpcWorkspace.state[215];
nmpcWorkspace.evGu[lRun1 * 56 + 6] = nmpcWorkspace.state[216];
nmpcWorkspace.evGu[lRun1 * 56 + 7] = nmpcWorkspace.state[217];
nmpcWorkspace.evGu[lRun1 * 56 + 8] = nmpcWorkspace.state[218];
nmpcWorkspace.evGu[lRun1 * 56 + 9] = nmpcWorkspace.state[219];
nmpcWorkspace.evGu[lRun1 * 56 + 10] = nmpcWorkspace.state[220];
nmpcWorkspace.evGu[lRun1 * 56 + 11] = nmpcWorkspace.state[221];
nmpcWorkspace.evGu[lRun1 * 56 + 12] = nmpcWorkspace.state[222];
nmpcWorkspace.evGu[lRun1 * 56 + 13] = nmpcWorkspace.state[223];
nmpcWorkspace.evGu[lRun1 * 56 + 14] = nmpcWorkspace.state[224];
nmpcWorkspace.evGu[lRun1 * 56 + 15] = nmpcWorkspace.state[225];
nmpcWorkspace.evGu[lRun1 * 56 + 16] = nmpcWorkspace.state[226];
nmpcWorkspace.evGu[lRun1 * 56 + 17] = nmpcWorkspace.state[227];
nmpcWorkspace.evGu[lRun1 * 56 + 18] = nmpcWorkspace.state[228];
nmpcWorkspace.evGu[lRun1 * 56 + 19] = nmpcWorkspace.state[229];
nmpcWorkspace.evGu[lRun1 * 56 + 20] = nmpcWorkspace.state[230];
nmpcWorkspace.evGu[lRun1 * 56 + 21] = nmpcWorkspace.state[231];
nmpcWorkspace.evGu[lRun1 * 56 + 22] = nmpcWorkspace.state[232];
nmpcWorkspace.evGu[lRun1 * 56 + 23] = nmpcWorkspace.state[233];
nmpcWorkspace.evGu[lRun1 * 56 + 24] = nmpcWorkspace.state[234];
nmpcWorkspace.evGu[lRun1 * 56 + 25] = nmpcWorkspace.state[235];
nmpcWorkspace.evGu[lRun1 * 56 + 26] = nmpcWorkspace.state[236];
nmpcWorkspace.evGu[lRun1 * 56 + 27] = nmpcWorkspace.state[237];
nmpcWorkspace.evGu[lRun1 * 56 + 28] = nmpcWorkspace.state[238];
nmpcWorkspace.evGu[lRun1 * 56 + 29] = nmpcWorkspace.state[239];
nmpcWorkspace.evGu[lRun1 * 56 + 30] = nmpcWorkspace.state[240];
nmpcWorkspace.evGu[lRun1 * 56 + 31] = nmpcWorkspace.state[241];
nmpcWorkspace.evGu[lRun1 * 56 + 32] = nmpcWorkspace.state[242];
nmpcWorkspace.evGu[lRun1 * 56 + 33] = nmpcWorkspace.state[243];
nmpcWorkspace.evGu[lRun1 * 56 + 34] = nmpcWorkspace.state[244];
nmpcWorkspace.evGu[lRun1 * 56 + 35] = nmpcWorkspace.state[245];
nmpcWorkspace.evGu[lRun1 * 56 + 36] = nmpcWorkspace.state[246];
nmpcWorkspace.evGu[lRun1 * 56 + 37] = nmpcWorkspace.state[247];
nmpcWorkspace.evGu[lRun1 * 56 + 38] = nmpcWorkspace.state[248];
nmpcWorkspace.evGu[lRun1 * 56 + 39] = nmpcWorkspace.state[249];
nmpcWorkspace.evGu[lRun1 * 56 + 40] = nmpcWorkspace.state[250];
nmpcWorkspace.evGu[lRun1 * 56 + 41] = nmpcWorkspace.state[251];
nmpcWorkspace.evGu[lRun1 * 56 + 42] = nmpcWorkspace.state[252];
nmpcWorkspace.evGu[lRun1 * 56 + 43] = nmpcWorkspace.state[253];
nmpcWorkspace.evGu[lRun1 * 56 + 44] = nmpcWorkspace.state[254];
nmpcWorkspace.evGu[lRun1 * 56 + 45] = nmpcWorkspace.state[255];
nmpcWorkspace.evGu[lRun1 * 56 + 46] = nmpcWorkspace.state[256];
nmpcWorkspace.evGu[lRun1 * 56 + 47] = nmpcWorkspace.state[257];
nmpcWorkspace.evGu[lRun1 * 56 + 48] = nmpcWorkspace.state[258];
nmpcWorkspace.evGu[lRun1 * 56 + 49] = nmpcWorkspace.state[259];
nmpcWorkspace.evGu[lRun1 * 56 + 50] = nmpcWorkspace.state[260];
nmpcWorkspace.evGu[lRun1 * 56 + 51] = nmpcWorkspace.state[261];
nmpcWorkspace.evGu[lRun1 * 56 + 52] = nmpcWorkspace.state[262];
nmpcWorkspace.evGu[lRun1 * 56 + 53] = nmpcWorkspace.state[263];
nmpcWorkspace.evGu[lRun1 * 56 + 54] = nmpcWorkspace.state[264];
nmpcWorkspace.evGu[lRun1 * 56 + 55] = nmpcWorkspace.state[265];
}
return ret;
}

void nmpc_evaluateLSQ(const real_t* in, real_t* out)
{
const real_t* xd = in;
const real_t* u = in + 14;
const real_t* od = in + 18;
/* Vector of auxiliary variables; number of elements: 6. */
real_t* a = nmpcWorkspace.objAuxVar;

/* Compute intermediate quantities: */
a[0] = (od[0]-od[3]);
a[1] = (xd[0]-od[3]);
a[2] = (od[1]-od[4]);
a[3] = (xd[0]-od[4]);
a[4] = ((a[0]*a[1])+(a[2]*a[3]));
a[5] = (a[0]+a[2]);

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
out[4] = xd[4];
out[5] = xd[5];
out[6] = xd[6];
out[7] = (a[4]-(real_t)(1.0000000000000000e+00));
out[8] = xd[7];
out[9] = u[0];
out[10] = u[1];
out[11] = u[2];
out[12] = u[3];
out[13] = (real_t)(1.0000000000000000e+00);
out[14] = (real_t)(0.0000000000000000e+00);
out[15] = (real_t)(0.0000000000000000e+00);
out[16] = (real_t)(0.0000000000000000e+00);
out[17] = (real_t)(0.0000000000000000e+00);
out[18] = (real_t)(0.0000000000000000e+00);
out[19] = (real_t)(0.0000000000000000e+00);
out[20] = (real_t)(0.0000000000000000e+00);
out[21] = (real_t)(0.0000000000000000e+00);
out[22] = (real_t)(0.0000000000000000e+00);
out[23] = (real_t)(0.0000000000000000e+00);
out[24] = (real_t)(0.0000000000000000e+00);
out[25] = (real_t)(0.0000000000000000e+00);
out[26] = (real_t)(0.0000000000000000e+00);
out[27] = (real_t)(0.0000000000000000e+00);
out[28] = (real_t)(1.0000000000000000e+00);
out[29] = (real_t)(0.0000000000000000e+00);
out[30] = (real_t)(0.0000000000000000e+00);
out[31] = (real_t)(0.0000000000000000e+00);
out[32] = (real_t)(0.0000000000000000e+00);
out[33] = (real_t)(0.0000000000000000e+00);
out[34] = (real_t)(0.0000000000000000e+00);
out[35] = (real_t)(0.0000000000000000e+00);
out[36] = (real_t)(0.0000000000000000e+00);
out[37] = (real_t)(0.0000000000000000e+00);
out[38] = (real_t)(0.0000000000000000e+00);
out[39] = (real_t)(0.0000000000000000e+00);
out[40] = (real_t)(0.0000000000000000e+00);
out[41] = (real_t)(0.0000000000000000e+00);
out[42] = (real_t)(0.0000000000000000e+00);
out[43] = (real_t)(1.0000000000000000e+00);
out[44] = (real_t)(0.0000000000000000e+00);
out[45] = (real_t)(0.0000000000000000e+00);
out[46] = (real_t)(0.0000000000000000e+00);
out[47] = (real_t)(0.0000000000000000e+00);
out[48] = (real_t)(0.0000000000000000e+00);
out[49] = (real_t)(0.0000000000000000e+00);
out[50] = (real_t)(0.0000000000000000e+00);
out[51] = (real_t)(0.0000000000000000e+00);
out[52] = (real_t)(0.0000000000000000e+00);
out[53] = (real_t)(0.0000000000000000e+00);
out[54] = (real_t)(0.0000000000000000e+00);
out[55] = (real_t)(0.0000000000000000e+00);
out[56] = (real_t)(0.0000000000000000e+00);
out[57] = (real_t)(0.0000000000000000e+00);
out[58] = (real_t)(1.0000000000000000e+00);
out[59] = (real_t)(0.0000000000000000e+00);
out[60] = (real_t)(0.0000000000000000e+00);
out[61] = (real_t)(0.0000000000000000e+00);
out[62] = (real_t)(0.0000000000000000e+00);
out[63] = (real_t)(0.0000000000000000e+00);
out[64] = (real_t)(0.0000000000000000e+00);
out[65] = (real_t)(0.0000000000000000e+00);
out[66] = (real_t)(0.0000000000000000e+00);
out[67] = (real_t)(0.0000000000000000e+00);
out[68] = (real_t)(0.0000000000000000e+00);
out[69] = (real_t)(0.0000000000000000e+00);
out[70] = (real_t)(0.0000000000000000e+00);
out[71] = (real_t)(0.0000000000000000e+00);
out[72] = (real_t)(0.0000000000000000e+00);
out[73] = (real_t)(1.0000000000000000e+00);
out[74] = (real_t)(0.0000000000000000e+00);
out[75] = (real_t)(0.0000000000000000e+00);
out[76] = (real_t)(0.0000000000000000e+00);
out[77] = (real_t)(0.0000000000000000e+00);
out[78] = (real_t)(0.0000000000000000e+00);
out[79] = (real_t)(0.0000000000000000e+00);
out[80] = (real_t)(0.0000000000000000e+00);
out[81] = (real_t)(0.0000000000000000e+00);
out[82] = (real_t)(0.0000000000000000e+00);
out[83] = (real_t)(0.0000000000000000e+00);
out[84] = (real_t)(0.0000000000000000e+00);
out[85] = (real_t)(0.0000000000000000e+00);
out[86] = (real_t)(0.0000000000000000e+00);
out[87] = (real_t)(0.0000000000000000e+00);
out[88] = (real_t)(1.0000000000000000e+00);
out[89] = (real_t)(0.0000000000000000e+00);
out[90] = (real_t)(0.0000000000000000e+00);
out[91] = (real_t)(0.0000000000000000e+00);
out[92] = (real_t)(0.0000000000000000e+00);
out[93] = (real_t)(0.0000000000000000e+00);
out[94] = (real_t)(0.0000000000000000e+00);
out[95] = (real_t)(0.0000000000000000e+00);
out[96] = (real_t)(0.0000000000000000e+00);
out[97] = (real_t)(0.0000000000000000e+00);
out[98] = (real_t)(0.0000000000000000e+00);
out[99] = (real_t)(0.0000000000000000e+00);
out[100] = (real_t)(0.0000000000000000e+00);
out[101] = (real_t)(0.0000000000000000e+00);
out[102] = (real_t)(0.0000000000000000e+00);
out[103] = (real_t)(1.0000000000000000e+00);
out[104] = (real_t)(0.0000000000000000e+00);
out[105] = (real_t)(0.0000000000000000e+00);
out[106] = (real_t)(0.0000000000000000e+00);
out[107] = (real_t)(0.0000000000000000e+00);
out[108] = (real_t)(0.0000000000000000e+00);
out[109] = (real_t)(0.0000000000000000e+00);
out[110] = (real_t)(0.0000000000000000e+00);
out[111] = a[5];
out[112] = (real_t)(0.0000000000000000e+00);
out[113] = (real_t)(0.0000000000000000e+00);
out[114] = (real_t)(0.0000000000000000e+00);
out[115] = (real_t)(0.0000000000000000e+00);
out[116] = (real_t)(0.0000000000000000e+00);
out[117] = (real_t)(0.0000000000000000e+00);
out[118] = (real_t)(0.0000000000000000e+00);
out[119] = (real_t)(0.0000000000000000e+00);
out[120] = (real_t)(0.0000000000000000e+00);
out[121] = (real_t)(0.0000000000000000e+00);
out[122] = (real_t)(0.0000000000000000e+00);
out[123] = (real_t)(0.0000000000000000e+00);
out[124] = (real_t)(0.0000000000000000e+00);
out[125] = (real_t)(0.0000000000000000e+00);
out[126] = (real_t)(0.0000000000000000e+00);
out[127] = (real_t)(0.0000000000000000e+00);
out[128] = (real_t)(0.0000000000000000e+00);
out[129] = (real_t)(0.0000000000000000e+00);
out[130] = (real_t)(0.0000000000000000e+00);
out[131] = (real_t)(0.0000000000000000e+00);
out[132] = (real_t)(1.0000000000000000e+00);
out[133] = (real_t)(0.0000000000000000e+00);
out[134] = (real_t)(0.0000000000000000e+00);
out[135] = (real_t)(0.0000000000000000e+00);
out[136] = (real_t)(0.0000000000000000e+00);
out[137] = (real_t)(0.0000000000000000e+00);
out[138] = (real_t)(0.0000000000000000e+00);
out[139] = (real_t)(0.0000000000000000e+00);
out[140] = (real_t)(0.0000000000000000e+00);
out[141] = (real_t)(0.0000000000000000e+00);
out[142] = (real_t)(0.0000000000000000e+00);
out[143] = (real_t)(0.0000000000000000e+00);
out[144] = (real_t)(0.0000000000000000e+00);
out[145] = (real_t)(0.0000000000000000e+00);
out[146] = (real_t)(0.0000000000000000e+00);
out[147] = (real_t)(0.0000000000000000e+00);
out[148] = (real_t)(0.0000000000000000e+00);
out[149] = (real_t)(0.0000000000000000e+00);
out[150] = (real_t)(0.0000000000000000e+00);
out[151] = (real_t)(0.0000000000000000e+00);
out[152] = (real_t)(0.0000000000000000e+00);
out[153] = (real_t)(0.0000000000000000e+00);
out[154] = (real_t)(0.0000000000000000e+00);
out[155] = (real_t)(0.0000000000000000e+00);
out[156] = (real_t)(0.0000000000000000e+00);
out[157] = (real_t)(0.0000000000000000e+00);
out[158] = (real_t)(0.0000000000000000e+00);
out[159] = (real_t)(0.0000000000000000e+00);
out[160] = (real_t)(0.0000000000000000e+00);
out[161] = (real_t)(0.0000000000000000e+00);
out[162] = (real_t)(0.0000000000000000e+00);
out[163] = (real_t)(0.0000000000000000e+00);
out[164] = (real_t)(0.0000000000000000e+00);
out[165] = (real_t)(0.0000000000000000e+00);
out[166] = (real_t)(0.0000000000000000e+00);
out[167] = (real_t)(0.0000000000000000e+00);
out[168] = (real_t)(0.0000000000000000e+00);
out[169] = (real_t)(0.0000000000000000e+00);
out[170] = (real_t)(0.0000000000000000e+00);
out[171] = (real_t)(0.0000000000000000e+00);
out[172] = (real_t)(0.0000000000000000e+00);
out[173] = (real_t)(0.0000000000000000e+00);
out[174] = (real_t)(0.0000000000000000e+00);
out[175] = (real_t)(0.0000000000000000e+00);
out[176] = (real_t)(0.0000000000000000e+00);
out[177] = (real_t)(0.0000000000000000e+00);
out[178] = (real_t)(0.0000000000000000e+00);
out[179] = (real_t)(0.0000000000000000e+00);
out[180] = (real_t)(0.0000000000000000e+00);
out[181] = (real_t)(0.0000000000000000e+00);
out[182] = (real_t)(0.0000000000000000e+00);
out[183] = (real_t)(0.0000000000000000e+00);
out[184] = (real_t)(0.0000000000000000e+00);
out[185] = (real_t)(0.0000000000000000e+00);
out[186] = (real_t)(0.0000000000000000e+00);
out[187] = (real_t)(0.0000000000000000e+00);
out[188] = (real_t)(0.0000000000000000e+00);
out[189] = (real_t)(0.0000000000000000e+00);
out[190] = (real_t)(0.0000000000000000e+00);
out[191] = (real_t)(0.0000000000000000e+00);
out[192] = (real_t)(0.0000000000000000e+00);
out[193] = (real_t)(0.0000000000000000e+00);
out[194] = (real_t)(0.0000000000000000e+00);
}

void nmpc_evaluateLSQEndTerm(const real_t* in, real_t* out)
{
const real_t* xd = in;

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
out[4] = xd[4];
out[5] = xd[5];
out[6] = xd[6];
out[7] = xd[7];
}

void nmpc_setObjQ1Q2( real_t* const tmpFx, real_t* const tmpObjS, real_t* const tmpQ1, real_t* const tmpQ2 )
{
tmpQ2[0] = + tmpFx[0]*tmpObjS[0] + tmpFx[14]*tmpObjS[13] + tmpFx[28]*tmpObjS[26] + tmpFx[42]*tmpObjS[39] + tmpFx[56]*tmpObjS[52] + tmpFx[70]*tmpObjS[65] + tmpFx[84]*tmpObjS[78] + tmpFx[98]*tmpObjS[91] + tmpFx[112]*tmpObjS[104] + tmpFx[126]*tmpObjS[117] + tmpFx[140]*tmpObjS[130] + tmpFx[154]*tmpObjS[143] + tmpFx[168]*tmpObjS[156];
tmpQ2[1] = + tmpFx[0]*tmpObjS[1] + tmpFx[14]*tmpObjS[14] + tmpFx[28]*tmpObjS[27] + tmpFx[42]*tmpObjS[40] + tmpFx[56]*tmpObjS[53] + tmpFx[70]*tmpObjS[66] + tmpFx[84]*tmpObjS[79] + tmpFx[98]*tmpObjS[92] + tmpFx[112]*tmpObjS[105] + tmpFx[126]*tmpObjS[118] + tmpFx[140]*tmpObjS[131] + tmpFx[154]*tmpObjS[144] + tmpFx[168]*tmpObjS[157];
tmpQ2[2] = + tmpFx[0]*tmpObjS[2] + tmpFx[14]*tmpObjS[15] + tmpFx[28]*tmpObjS[28] + tmpFx[42]*tmpObjS[41] + tmpFx[56]*tmpObjS[54] + tmpFx[70]*tmpObjS[67] + tmpFx[84]*tmpObjS[80] + tmpFx[98]*tmpObjS[93] + tmpFx[112]*tmpObjS[106] + tmpFx[126]*tmpObjS[119] + tmpFx[140]*tmpObjS[132] + tmpFx[154]*tmpObjS[145] + tmpFx[168]*tmpObjS[158];
tmpQ2[3] = + tmpFx[0]*tmpObjS[3] + tmpFx[14]*tmpObjS[16] + tmpFx[28]*tmpObjS[29] + tmpFx[42]*tmpObjS[42] + tmpFx[56]*tmpObjS[55] + tmpFx[70]*tmpObjS[68] + tmpFx[84]*tmpObjS[81] + tmpFx[98]*tmpObjS[94] + tmpFx[112]*tmpObjS[107] + tmpFx[126]*tmpObjS[120] + tmpFx[140]*tmpObjS[133] + tmpFx[154]*tmpObjS[146] + tmpFx[168]*tmpObjS[159];
tmpQ2[4] = + tmpFx[0]*tmpObjS[4] + tmpFx[14]*tmpObjS[17] + tmpFx[28]*tmpObjS[30] + tmpFx[42]*tmpObjS[43] + tmpFx[56]*tmpObjS[56] + tmpFx[70]*tmpObjS[69] + tmpFx[84]*tmpObjS[82] + tmpFx[98]*tmpObjS[95] + tmpFx[112]*tmpObjS[108] + tmpFx[126]*tmpObjS[121] + tmpFx[140]*tmpObjS[134] + tmpFx[154]*tmpObjS[147] + tmpFx[168]*tmpObjS[160];
tmpQ2[5] = + tmpFx[0]*tmpObjS[5] + tmpFx[14]*tmpObjS[18] + tmpFx[28]*tmpObjS[31] + tmpFx[42]*tmpObjS[44] + tmpFx[56]*tmpObjS[57] + tmpFx[70]*tmpObjS[70] + tmpFx[84]*tmpObjS[83] + tmpFx[98]*tmpObjS[96] + tmpFx[112]*tmpObjS[109] + tmpFx[126]*tmpObjS[122] + tmpFx[140]*tmpObjS[135] + tmpFx[154]*tmpObjS[148] + tmpFx[168]*tmpObjS[161];
tmpQ2[6] = + tmpFx[0]*tmpObjS[6] + tmpFx[14]*tmpObjS[19] + tmpFx[28]*tmpObjS[32] + tmpFx[42]*tmpObjS[45] + tmpFx[56]*tmpObjS[58] + tmpFx[70]*tmpObjS[71] + tmpFx[84]*tmpObjS[84] + tmpFx[98]*tmpObjS[97] + tmpFx[112]*tmpObjS[110] + tmpFx[126]*tmpObjS[123] + tmpFx[140]*tmpObjS[136] + tmpFx[154]*tmpObjS[149] + tmpFx[168]*tmpObjS[162];
tmpQ2[7] = + tmpFx[0]*tmpObjS[7] + tmpFx[14]*tmpObjS[20] + tmpFx[28]*tmpObjS[33] + tmpFx[42]*tmpObjS[46] + tmpFx[56]*tmpObjS[59] + tmpFx[70]*tmpObjS[72] + tmpFx[84]*tmpObjS[85] + tmpFx[98]*tmpObjS[98] + tmpFx[112]*tmpObjS[111] + tmpFx[126]*tmpObjS[124] + tmpFx[140]*tmpObjS[137] + tmpFx[154]*tmpObjS[150] + tmpFx[168]*tmpObjS[163];
tmpQ2[8] = + tmpFx[0]*tmpObjS[8] + tmpFx[14]*tmpObjS[21] + tmpFx[28]*tmpObjS[34] + tmpFx[42]*tmpObjS[47] + tmpFx[56]*tmpObjS[60] + tmpFx[70]*tmpObjS[73] + tmpFx[84]*tmpObjS[86] + tmpFx[98]*tmpObjS[99] + tmpFx[112]*tmpObjS[112] + tmpFx[126]*tmpObjS[125] + tmpFx[140]*tmpObjS[138] + tmpFx[154]*tmpObjS[151] + tmpFx[168]*tmpObjS[164];
tmpQ2[9] = + tmpFx[0]*tmpObjS[9] + tmpFx[14]*tmpObjS[22] + tmpFx[28]*tmpObjS[35] + tmpFx[42]*tmpObjS[48] + tmpFx[56]*tmpObjS[61] + tmpFx[70]*tmpObjS[74] + tmpFx[84]*tmpObjS[87] + tmpFx[98]*tmpObjS[100] + tmpFx[112]*tmpObjS[113] + tmpFx[126]*tmpObjS[126] + tmpFx[140]*tmpObjS[139] + tmpFx[154]*tmpObjS[152] + tmpFx[168]*tmpObjS[165];
tmpQ2[10] = + tmpFx[0]*tmpObjS[10] + tmpFx[14]*tmpObjS[23] + tmpFx[28]*tmpObjS[36] + tmpFx[42]*tmpObjS[49] + tmpFx[56]*tmpObjS[62] + tmpFx[70]*tmpObjS[75] + tmpFx[84]*tmpObjS[88] + tmpFx[98]*tmpObjS[101] + tmpFx[112]*tmpObjS[114] + tmpFx[126]*tmpObjS[127] + tmpFx[140]*tmpObjS[140] + tmpFx[154]*tmpObjS[153] + tmpFx[168]*tmpObjS[166];
tmpQ2[11] = + tmpFx[0]*tmpObjS[11] + tmpFx[14]*tmpObjS[24] + tmpFx[28]*tmpObjS[37] + tmpFx[42]*tmpObjS[50] + tmpFx[56]*tmpObjS[63] + tmpFx[70]*tmpObjS[76] + tmpFx[84]*tmpObjS[89] + tmpFx[98]*tmpObjS[102] + tmpFx[112]*tmpObjS[115] + tmpFx[126]*tmpObjS[128] + tmpFx[140]*tmpObjS[141] + tmpFx[154]*tmpObjS[154] + tmpFx[168]*tmpObjS[167];
tmpQ2[12] = + tmpFx[0]*tmpObjS[12] + tmpFx[14]*tmpObjS[25] + tmpFx[28]*tmpObjS[38] + tmpFx[42]*tmpObjS[51] + tmpFx[56]*tmpObjS[64] + tmpFx[70]*tmpObjS[77] + tmpFx[84]*tmpObjS[90] + tmpFx[98]*tmpObjS[103] + tmpFx[112]*tmpObjS[116] + tmpFx[126]*tmpObjS[129] + tmpFx[140]*tmpObjS[142] + tmpFx[154]*tmpObjS[155] + tmpFx[168]*tmpObjS[168];
tmpQ2[13] = + tmpFx[1]*tmpObjS[0] + tmpFx[15]*tmpObjS[13] + tmpFx[29]*tmpObjS[26] + tmpFx[43]*tmpObjS[39] + tmpFx[57]*tmpObjS[52] + tmpFx[71]*tmpObjS[65] + tmpFx[85]*tmpObjS[78] + tmpFx[99]*tmpObjS[91] + tmpFx[113]*tmpObjS[104] + tmpFx[127]*tmpObjS[117] + tmpFx[141]*tmpObjS[130] + tmpFx[155]*tmpObjS[143] + tmpFx[169]*tmpObjS[156];
tmpQ2[14] = + tmpFx[1]*tmpObjS[1] + tmpFx[15]*tmpObjS[14] + tmpFx[29]*tmpObjS[27] + tmpFx[43]*tmpObjS[40] + tmpFx[57]*tmpObjS[53] + tmpFx[71]*tmpObjS[66] + tmpFx[85]*tmpObjS[79] + tmpFx[99]*tmpObjS[92] + tmpFx[113]*tmpObjS[105] + tmpFx[127]*tmpObjS[118] + tmpFx[141]*tmpObjS[131] + tmpFx[155]*tmpObjS[144] + tmpFx[169]*tmpObjS[157];
tmpQ2[15] = + tmpFx[1]*tmpObjS[2] + tmpFx[15]*tmpObjS[15] + tmpFx[29]*tmpObjS[28] + tmpFx[43]*tmpObjS[41] + tmpFx[57]*tmpObjS[54] + tmpFx[71]*tmpObjS[67] + tmpFx[85]*tmpObjS[80] + tmpFx[99]*tmpObjS[93] + tmpFx[113]*tmpObjS[106] + tmpFx[127]*tmpObjS[119] + tmpFx[141]*tmpObjS[132] + tmpFx[155]*tmpObjS[145] + tmpFx[169]*tmpObjS[158];
tmpQ2[16] = + tmpFx[1]*tmpObjS[3] + tmpFx[15]*tmpObjS[16] + tmpFx[29]*tmpObjS[29] + tmpFx[43]*tmpObjS[42] + tmpFx[57]*tmpObjS[55] + tmpFx[71]*tmpObjS[68] + tmpFx[85]*tmpObjS[81] + tmpFx[99]*tmpObjS[94] + tmpFx[113]*tmpObjS[107] + tmpFx[127]*tmpObjS[120] + tmpFx[141]*tmpObjS[133] + tmpFx[155]*tmpObjS[146] + tmpFx[169]*tmpObjS[159];
tmpQ2[17] = + tmpFx[1]*tmpObjS[4] + tmpFx[15]*tmpObjS[17] + tmpFx[29]*tmpObjS[30] + tmpFx[43]*tmpObjS[43] + tmpFx[57]*tmpObjS[56] + tmpFx[71]*tmpObjS[69] + tmpFx[85]*tmpObjS[82] + tmpFx[99]*tmpObjS[95] + tmpFx[113]*tmpObjS[108] + tmpFx[127]*tmpObjS[121] + tmpFx[141]*tmpObjS[134] + tmpFx[155]*tmpObjS[147] + tmpFx[169]*tmpObjS[160];
tmpQ2[18] = + tmpFx[1]*tmpObjS[5] + tmpFx[15]*tmpObjS[18] + tmpFx[29]*tmpObjS[31] + tmpFx[43]*tmpObjS[44] + tmpFx[57]*tmpObjS[57] + tmpFx[71]*tmpObjS[70] + tmpFx[85]*tmpObjS[83] + tmpFx[99]*tmpObjS[96] + tmpFx[113]*tmpObjS[109] + tmpFx[127]*tmpObjS[122] + tmpFx[141]*tmpObjS[135] + tmpFx[155]*tmpObjS[148] + tmpFx[169]*tmpObjS[161];
tmpQ2[19] = + tmpFx[1]*tmpObjS[6] + tmpFx[15]*tmpObjS[19] + tmpFx[29]*tmpObjS[32] + tmpFx[43]*tmpObjS[45] + tmpFx[57]*tmpObjS[58] + tmpFx[71]*tmpObjS[71] + tmpFx[85]*tmpObjS[84] + tmpFx[99]*tmpObjS[97] + tmpFx[113]*tmpObjS[110] + tmpFx[127]*tmpObjS[123] + tmpFx[141]*tmpObjS[136] + tmpFx[155]*tmpObjS[149] + tmpFx[169]*tmpObjS[162];
tmpQ2[20] = + tmpFx[1]*tmpObjS[7] + tmpFx[15]*tmpObjS[20] + tmpFx[29]*tmpObjS[33] + tmpFx[43]*tmpObjS[46] + tmpFx[57]*tmpObjS[59] + tmpFx[71]*tmpObjS[72] + tmpFx[85]*tmpObjS[85] + tmpFx[99]*tmpObjS[98] + tmpFx[113]*tmpObjS[111] + tmpFx[127]*tmpObjS[124] + tmpFx[141]*tmpObjS[137] + tmpFx[155]*tmpObjS[150] + tmpFx[169]*tmpObjS[163];
tmpQ2[21] = + tmpFx[1]*tmpObjS[8] + tmpFx[15]*tmpObjS[21] + tmpFx[29]*tmpObjS[34] + tmpFx[43]*tmpObjS[47] + tmpFx[57]*tmpObjS[60] + tmpFx[71]*tmpObjS[73] + tmpFx[85]*tmpObjS[86] + tmpFx[99]*tmpObjS[99] + tmpFx[113]*tmpObjS[112] + tmpFx[127]*tmpObjS[125] + tmpFx[141]*tmpObjS[138] + tmpFx[155]*tmpObjS[151] + tmpFx[169]*tmpObjS[164];
tmpQ2[22] = + tmpFx[1]*tmpObjS[9] + tmpFx[15]*tmpObjS[22] + tmpFx[29]*tmpObjS[35] + tmpFx[43]*tmpObjS[48] + tmpFx[57]*tmpObjS[61] + tmpFx[71]*tmpObjS[74] + tmpFx[85]*tmpObjS[87] + tmpFx[99]*tmpObjS[100] + tmpFx[113]*tmpObjS[113] + tmpFx[127]*tmpObjS[126] + tmpFx[141]*tmpObjS[139] + tmpFx[155]*tmpObjS[152] + tmpFx[169]*tmpObjS[165];
tmpQ2[23] = + tmpFx[1]*tmpObjS[10] + tmpFx[15]*tmpObjS[23] + tmpFx[29]*tmpObjS[36] + tmpFx[43]*tmpObjS[49] + tmpFx[57]*tmpObjS[62] + tmpFx[71]*tmpObjS[75] + tmpFx[85]*tmpObjS[88] + tmpFx[99]*tmpObjS[101] + tmpFx[113]*tmpObjS[114] + tmpFx[127]*tmpObjS[127] + tmpFx[141]*tmpObjS[140] + tmpFx[155]*tmpObjS[153] + tmpFx[169]*tmpObjS[166];
tmpQ2[24] = + tmpFx[1]*tmpObjS[11] + tmpFx[15]*tmpObjS[24] + tmpFx[29]*tmpObjS[37] + tmpFx[43]*tmpObjS[50] + tmpFx[57]*tmpObjS[63] + tmpFx[71]*tmpObjS[76] + tmpFx[85]*tmpObjS[89] + tmpFx[99]*tmpObjS[102] + tmpFx[113]*tmpObjS[115] + tmpFx[127]*tmpObjS[128] + tmpFx[141]*tmpObjS[141] + tmpFx[155]*tmpObjS[154] + tmpFx[169]*tmpObjS[167];
tmpQ2[25] = + tmpFx[1]*tmpObjS[12] + tmpFx[15]*tmpObjS[25] + tmpFx[29]*tmpObjS[38] + tmpFx[43]*tmpObjS[51] + tmpFx[57]*tmpObjS[64] + tmpFx[71]*tmpObjS[77] + tmpFx[85]*tmpObjS[90] + tmpFx[99]*tmpObjS[103] + tmpFx[113]*tmpObjS[116] + tmpFx[127]*tmpObjS[129] + tmpFx[141]*tmpObjS[142] + tmpFx[155]*tmpObjS[155] + tmpFx[169]*tmpObjS[168];
tmpQ2[26] = + tmpFx[2]*tmpObjS[0] + tmpFx[16]*tmpObjS[13] + tmpFx[30]*tmpObjS[26] + tmpFx[44]*tmpObjS[39] + tmpFx[58]*tmpObjS[52] + tmpFx[72]*tmpObjS[65] + tmpFx[86]*tmpObjS[78] + tmpFx[100]*tmpObjS[91] + tmpFx[114]*tmpObjS[104] + tmpFx[128]*tmpObjS[117] + tmpFx[142]*tmpObjS[130] + tmpFx[156]*tmpObjS[143] + tmpFx[170]*tmpObjS[156];
tmpQ2[27] = + tmpFx[2]*tmpObjS[1] + tmpFx[16]*tmpObjS[14] + tmpFx[30]*tmpObjS[27] + tmpFx[44]*tmpObjS[40] + tmpFx[58]*tmpObjS[53] + tmpFx[72]*tmpObjS[66] + tmpFx[86]*tmpObjS[79] + tmpFx[100]*tmpObjS[92] + tmpFx[114]*tmpObjS[105] + tmpFx[128]*tmpObjS[118] + tmpFx[142]*tmpObjS[131] + tmpFx[156]*tmpObjS[144] + tmpFx[170]*tmpObjS[157];
tmpQ2[28] = + tmpFx[2]*tmpObjS[2] + tmpFx[16]*tmpObjS[15] + tmpFx[30]*tmpObjS[28] + tmpFx[44]*tmpObjS[41] + tmpFx[58]*tmpObjS[54] + tmpFx[72]*tmpObjS[67] + tmpFx[86]*tmpObjS[80] + tmpFx[100]*tmpObjS[93] + tmpFx[114]*tmpObjS[106] + tmpFx[128]*tmpObjS[119] + tmpFx[142]*tmpObjS[132] + tmpFx[156]*tmpObjS[145] + tmpFx[170]*tmpObjS[158];
tmpQ2[29] = + tmpFx[2]*tmpObjS[3] + tmpFx[16]*tmpObjS[16] + tmpFx[30]*tmpObjS[29] + tmpFx[44]*tmpObjS[42] + tmpFx[58]*tmpObjS[55] + tmpFx[72]*tmpObjS[68] + tmpFx[86]*tmpObjS[81] + tmpFx[100]*tmpObjS[94] + tmpFx[114]*tmpObjS[107] + tmpFx[128]*tmpObjS[120] + tmpFx[142]*tmpObjS[133] + tmpFx[156]*tmpObjS[146] + tmpFx[170]*tmpObjS[159];
tmpQ2[30] = + tmpFx[2]*tmpObjS[4] + tmpFx[16]*tmpObjS[17] + tmpFx[30]*tmpObjS[30] + tmpFx[44]*tmpObjS[43] + tmpFx[58]*tmpObjS[56] + tmpFx[72]*tmpObjS[69] + tmpFx[86]*tmpObjS[82] + tmpFx[100]*tmpObjS[95] + tmpFx[114]*tmpObjS[108] + tmpFx[128]*tmpObjS[121] + tmpFx[142]*tmpObjS[134] + tmpFx[156]*tmpObjS[147] + tmpFx[170]*tmpObjS[160];
tmpQ2[31] = + tmpFx[2]*tmpObjS[5] + tmpFx[16]*tmpObjS[18] + tmpFx[30]*tmpObjS[31] + tmpFx[44]*tmpObjS[44] + tmpFx[58]*tmpObjS[57] + tmpFx[72]*tmpObjS[70] + tmpFx[86]*tmpObjS[83] + tmpFx[100]*tmpObjS[96] + tmpFx[114]*tmpObjS[109] + tmpFx[128]*tmpObjS[122] + tmpFx[142]*tmpObjS[135] + tmpFx[156]*tmpObjS[148] + tmpFx[170]*tmpObjS[161];
tmpQ2[32] = + tmpFx[2]*tmpObjS[6] + tmpFx[16]*tmpObjS[19] + tmpFx[30]*tmpObjS[32] + tmpFx[44]*tmpObjS[45] + tmpFx[58]*tmpObjS[58] + tmpFx[72]*tmpObjS[71] + tmpFx[86]*tmpObjS[84] + tmpFx[100]*tmpObjS[97] + tmpFx[114]*tmpObjS[110] + tmpFx[128]*tmpObjS[123] + tmpFx[142]*tmpObjS[136] + tmpFx[156]*tmpObjS[149] + tmpFx[170]*tmpObjS[162];
tmpQ2[33] = + tmpFx[2]*tmpObjS[7] + tmpFx[16]*tmpObjS[20] + tmpFx[30]*tmpObjS[33] + tmpFx[44]*tmpObjS[46] + tmpFx[58]*tmpObjS[59] + tmpFx[72]*tmpObjS[72] + tmpFx[86]*tmpObjS[85] + tmpFx[100]*tmpObjS[98] + tmpFx[114]*tmpObjS[111] + tmpFx[128]*tmpObjS[124] + tmpFx[142]*tmpObjS[137] + tmpFx[156]*tmpObjS[150] + tmpFx[170]*tmpObjS[163];
tmpQ2[34] = + tmpFx[2]*tmpObjS[8] + tmpFx[16]*tmpObjS[21] + tmpFx[30]*tmpObjS[34] + tmpFx[44]*tmpObjS[47] + tmpFx[58]*tmpObjS[60] + tmpFx[72]*tmpObjS[73] + tmpFx[86]*tmpObjS[86] + tmpFx[100]*tmpObjS[99] + tmpFx[114]*tmpObjS[112] + tmpFx[128]*tmpObjS[125] + tmpFx[142]*tmpObjS[138] + tmpFx[156]*tmpObjS[151] + tmpFx[170]*tmpObjS[164];
tmpQ2[35] = + tmpFx[2]*tmpObjS[9] + tmpFx[16]*tmpObjS[22] + tmpFx[30]*tmpObjS[35] + tmpFx[44]*tmpObjS[48] + tmpFx[58]*tmpObjS[61] + tmpFx[72]*tmpObjS[74] + tmpFx[86]*tmpObjS[87] + tmpFx[100]*tmpObjS[100] + tmpFx[114]*tmpObjS[113] + tmpFx[128]*tmpObjS[126] + tmpFx[142]*tmpObjS[139] + tmpFx[156]*tmpObjS[152] + tmpFx[170]*tmpObjS[165];
tmpQ2[36] = + tmpFx[2]*tmpObjS[10] + tmpFx[16]*tmpObjS[23] + tmpFx[30]*tmpObjS[36] + tmpFx[44]*tmpObjS[49] + tmpFx[58]*tmpObjS[62] + tmpFx[72]*tmpObjS[75] + tmpFx[86]*tmpObjS[88] + tmpFx[100]*tmpObjS[101] + tmpFx[114]*tmpObjS[114] + tmpFx[128]*tmpObjS[127] + tmpFx[142]*tmpObjS[140] + tmpFx[156]*tmpObjS[153] + tmpFx[170]*tmpObjS[166];
tmpQ2[37] = + tmpFx[2]*tmpObjS[11] + tmpFx[16]*tmpObjS[24] + tmpFx[30]*tmpObjS[37] + tmpFx[44]*tmpObjS[50] + tmpFx[58]*tmpObjS[63] + tmpFx[72]*tmpObjS[76] + tmpFx[86]*tmpObjS[89] + tmpFx[100]*tmpObjS[102] + tmpFx[114]*tmpObjS[115] + tmpFx[128]*tmpObjS[128] + tmpFx[142]*tmpObjS[141] + tmpFx[156]*tmpObjS[154] + tmpFx[170]*tmpObjS[167];
tmpQ2[38] = + tmpFx[2]*tmpObjS[12] + tmpFx[16]*tmpObjS[25] + tmpFx[30]*tmpObjS[38] + tmpFx[44]*tmpObjS[51] + tmpFx[58]*tmpObjS[64] + tmpFx[72]*tmpObjS[77] + tmpFx[86]*tmpObjS[90] + tmpFx[100]*tmpObjS[103] + tmpFx[114]*tmpObjS[116] + tmpFx[128]*tmpObjS[129] + tmpFx[142]*tmpObjS[142] + tmpFx[156]*tmpObjS[155] + tmpFx[170]*tmpObjS[168];
tmpQ2[39] = + tmpFx[3]*tmpObjS[0] + tmpFx[17]*tmpObjS[13] + tmpFx[31]*tmpObjS[26] + tmpFx[45]*tmpObjS[39] + tmpFx[59]*tmpObjS[52] + tmpFx[73]*tmpObjS[65] + tmpFx[87]*tmpObjS[78] + tmpFx[101]*tmpObjS[91] + tmpFx[115]*tmpObjS[104] + tmpFx[129]*tmpObjS[117] + tmpFx[143]*tmpObjS[130] + tmpFx[157]*tmpObjS[143] + tmpFx[171]*tmpObjS[156];
tmpQ2[40] = + tmpFx[3]*tmpObjS[1] + tmpFx[17]*tmpObjS[14] + tmpFx[31]*tmpObjS[27] + tmpFx[45]*tmpObjS[40] + tmpFx[59]*tmpObjS[53] + tmpFx[73]*tmpObjS[66] + tmpFx[87]*tmpObjS[79] + tmpFx[101]*tmpObjS[92] + tmpFx[115]*tmpObjS[105] + tmpFx[129]*tmpObjS[118] + tmpFx[143]*tmpObjS[131] + tmpFx[157]*tmpObjS[144] + tmpFx[171]*tmpObjS[157];
tmpQ2[41] = + tmpFx[3]*tmpObjS[2] + tmpFx[17]*tmpObjS[15] + tmpFx[31]*tmpObjS[28] + tmpFx[45]*tmpObjS[41] + tmpFx[59]*tmpObjS[54] + tmpFx[73]*tmpObjS[67] + tmpFx[87]*tmpObjS[80] + tmpFx[101]*tmpObjS[93] + tmpFx[115]*tmpObjS[106] + tmpFx[129]*tmpObjS[119] + tmpFx[143]*tmpObjS[132] + tmpFx[157]*tmpObjS[145] + tmpFx[171]*tmpObjS[158];
tmpQ2[42] = + tmpFx[3]*tmpObjS[3] + tmpFx[17]*tmpObjS[16] + tmpFx[31]*tmpObjS[29] + tmpFx[45]*tmpObjS[42] + tmpFx[59]*tmpObjS[55] + tmpFx[73]*tmpObjS[68] + tmpFx[87]*tmpObjS[81] + tmpFx[101]*tmpObjS[94] + tmpFx[115]*tmpObjS[107] + tmpFx[129]*tmpObjS[120] + tmpFx[143]*tmpObjS[133] + tmpFx[157]*tmpObjS[146] + tmpFx[171]*tmpObjS[159];
tmpQ2[43] = + tmpFx[3]*tmpObjS[4] + tmpFx[17]*tmpObjS[17] + tmpFx[31]*tmpObjS[30] + tmpFx[45]*tmpObjS[43] + tmpFx[59]*tmpObjS[56] + tmpFx[73]*tmpObjS[69] + tmpFx[87]*tmpObjS[82] + tmpFx[101]*tmpObjS[95] + tmpFx[115]*tmpObjS[108] + tmpFx[129]*tmpObjS[121] + tmpFx[143]*tmpObjS[134] + tmpFx[157]*tmpObjS[147] + tmpFx[171]*tmpObjS[160];
tmpQ2[44] = + tmpFx[3]*tmpObjS[5] + tmpFx[17]*tmpObjS[18] + tmpFx[31]*tmpObjS[31] + tmpFx[45]*tmpObjS[44] + tmpFx[59]*tmpObjS[57] + tmpFx[73]*tmpObjS[70] + tmpFx[87]*tmpObjS[83] + tmpFx[101]*tmpObjS[96] + tmpFx[115]*tmpObjS[109] + tmpFx[129]*tmpObjS[122] + tmpFx[143]*tmpObjS[135] + tmpFx[157]*tmpObjS[148] + tmpFx[171]*tmpObjS[161];
tmpQ2[45] = + tmpFx[3]*tmpObjS[6] + tmpFx[17]*tmpObjS[19] + tmpFx[31]*tmpObjS[32] + tmpFx[45]*tmpObjS[45] + tmpFx[59]*tmpObjS[58] + tmpFx[73]*tmpObjS[71] + tmpFx[87]*tmpObjS[84] + tmpFx[101]*tmpObjS[97] + tmpFx[115]*tmpObjS[110] + tmpFx[129]*tmpObjS[123] + tmpFx[143]*tmpObjS[136] + tmpFx[157]*tmpObjS[149] + tmpFx[171]*tmpObjS[162];
tmpQ2[46] = + tmpFx[3]*tmpObjS[7] + tmpFx[17]*tmpObjS[20] + tmpFx[31]*tmpObjS[33] + tmpFx[45]*tmpObjS[46] + tmpFx[59]*tmpObjS[59] + tmpFx[73]*tmpObjS[72] + tmpFx[87]*tmpObjS[85] + tmpFx[101]*tmpObjS[98] + tmpFx[115]*tmpObjS[111] + tmpFx[129]*tmpObjS[124] + tmpFx[143]*tmpObjS[137] + tmpFx[157]*tmpObjS[150] + tmpFx[171]*tmpObjS[163];
tmpQ2[47] = + tmpFx[3]*tmpObjS[8] + tmpFx[17]*tmpObjS[21] + tmpFx[31]*tmpObjS[34] + tmpFx[45]*tmpObjS[47] + tmpFx[59]*tmpObjS[60] + tmpFx[73]*tmpObjS[73] + tmpFx[87]*tmpObjS[86] + tmpFx[101]*tmpObjS[99] + tmpFx[115]*tmpObjS[112] + tmpFx[129]*tmpObjS[125] + tmpFx[143]*tmpObjS[138] + tmpFx[157]*tmpObjS[151] + tmpFx[171]*tmpObjS[164];
tmpQ2[48] = + tmpFx[3]*tmpObjS[9] + tmpFx[17]*tmpObjS[22] + tmpFx[31]*tmpObjS[35] + tmpFx[45]*tmpObjS[48] + tmpFx[59]*tmpObjS[61] + tmpFx[73]*tmpObjS[74] + tmpFx[87]*tmpObjS[87] + tmpFx[101]*tmpObjS[100] + tmpFx[115]*tmpObjS[113] + tmpFx[129]*tmpObjS[126] + tmpFx[143]*tmpObjS[139] + tmpFx[157]*tmpObjS[152] + tmpFx[171]*tmpObjS[165];
tmpQ2[49] = + tmpFx[3]*tmpObjS[10] + tmpFx[17]*tmpObjS[23] + tmpFx[31]*tmpObjS[36] + tmpFx[45]*tmpObjS[49] + tmpFx[59]*tmpObjS[62] + tmpFx[73]*tmpObjS[75] + tmpFx[87]*tmpObjS[88] + tmpFx[101]*tmpObjS[101] + tmpFx[115]*tmpObjS[114] + tmpFx[129]*tmpObjS[127] + tmpFx[143]*tmpObjS[140] + tmpFx[157]*tmpObjS[153] + tmpFx[171]*tmpObjS[166];
tmpQ2[50] = + tmpFx[3]*tmpObjS[11] + tmpFx[17]*tmpObjS[24] + tmpFx[31]*tmpObjS[37] + tmpFx[45]*tmpObjS[50] + tmpFx[59]*tmpObjS[63] + tmpFx[73]*tmpObjS[76] + tmpFx[87]*tmpObjS[89] + tmpFx[101]*tmpObjS[102] + tmpFx[115]*tmpObjS[115] + tmpFx[129]*tmpObjS[128] + tmpFx[143]*tmpObjS[141] + tmpFx[157]*tmpObjS[154] + tmpFx[171]*tmpObjS[167];
tmpQ2[51] = + tmpFx[3]*tmpObjS[12] + tmpFx[17]*tmpObjS[25] + tmpFx[31]*tmpObjS[38] + tmpFx[45]*tmpObjS[51] + tmpFx[59]*tmpObjS[64] + tmpFx[73]*tmpObjS[77] + tmpFx[87]*tmpObjS[90] + tmpFx[101]*tmpObjS[103] + tmpFx[115]*tmpObjS[116] + tmpFx[129]*tmpObjS[129] + tmpFx[143]*tmpObjS[142] + tmpFx[157]*tmpObjS[155] + tmpFx[171]*tmpObjS[168];
tmpQ2[52] = + tmpFx[4]*tmpObjS[0] + tmpFx[18]*tmpObjS[13] + tmpFx[32]*tmpObjS[26] + tmpFx[46]*tmpObjS[39] + tmpFx[60]*tmpObjS[52] + tmpFx[74]*tmpObjS[65] + tmpFx[88]*tmpObjS[78] + tmpFx[102]*tmpObjS[91] + tmpFx[116]*tmpObjS[104] + tmpFx[130]*tmpObjS[117] + tmpFx[144]*tmpObjS[130] + tmpFx[158]*tmpObjS[143] + tmpFx[172]*tmpObjS[156];
tmpQ2[53] = + tmpFx[4]*tmpObjS[1] + tmpFx[18]*tmpObjS[14] + tmpFx[32]*tmpObjS[27] + tmpFx[46]*tmpObjS[40] + tmpFx[60]*tmpObjS[53] + tmpFx[74]*tmpObjS[66] + tmpFx[88]*tmpObjS[79] + tmpFx[102]*tmpObjS[92] + tmpFx[116]*tmpObjS[105] + tmpFx[130]*tmpObjS[118] + tmpFx[144]*tmpObjS[131] + tmpFx[158]*tmpObjS[144] + tmpFx[172]*tmpObjS[157];
tmpQ2[54] = + tmpFx[4]*tmpObjS[2] + tmpFx[18]*tmpObjS[15] + tmpFx[32]*tmpObjS[28] + tmpFx[46]*tmpObjS[41] + tmpFx[60]*tmpObjS[54] + tmpFx[74]*tmpObjS[67] + tmpFx[88]*tmpObjS[80] + tmpFx[102]*tmpObjS[93] + tmpFx[116]*tmpObjS[106] + tmpFx[130]*tmpObjS[119] + tmpFx[144]*tmpObjS[132] + tmpFx[158]*tmpObjS[145] + tmpFx[172]*tmpObjS[158];
tmpQ2[55] = + tmpFx[4]*tmpObjS[3] + tmpFx[18]*tmpObjS[16] + tmpFx[32]*tmpObjS[29] + tmpFx[46]*tmpObjS[42] + tmpFx[60]*tmpObjS[55] + tmpFx[74]*tmpObjS[68] + tmpFx[88]*tmpObjS[81] + tmpFx[102]*tmpObjS[94] + tmpFx[116]*tmpObjS[107] + tmpFx[130]*tmpObjS[120] + tmpFx[144]*tmpObjS[133] + tmpFx[158]*tmpObjS[146] + tmpFx[172]*tmpObjS[159];
tmpQ2[56] = + tmpFx[4]*tmpObjS[4] + tmpFx[18]*tmpObjS[17] + tmpFx[32]*tmpObjS[30] + tmpFx[46]*tmpObjS[43] + tmpFx[60]*tmpObjS[56] + tmpFx[74]*tmpObjS[69] + tmpFx[88]*tmpObjS[82] + tmpFx[102]*tmpObjS[95] + tmpFx[116]*tmpObjS[108] + tmpFx[130]*tmpObjS[121] + tmpFx[144]*tmpObjS[134] + tmpFx[158]*tmpObjS[147] + tmpFx[172]*tmpObjS[160];
tmpQ2[57] = + tmpFx[4]*tmpObjS[5] + tmpFx[18]*tmpObjS[18] + tmpFx[32]*tmpObjS[31] + tmpFx[46]*tmpObjS[44] + tmpFx[60]*tmpObjS[57] + tmpFx[74]*tmpObjS[70] + tmpFx[88]*tmpObjS[83] + tmpFx[102]*tmpObjS[96] + tmpFx[116]*tmpObjS[109] + tmpFx[130]*tmpObjS[122] + tmpFx[144]*tmpObjS[135] + tmpFx[158]*tmpObjS[148] + tmpFx[172]*tmpObjS[161];
tmpQ2[58] = + tmpFx[4]*tmpObjS[6] + tmpFx[18]*tmpObjS[19] + tmpFx[32]*tmpObjS[32] + tmpFx[46]*tmpObjS[45] + tmpFx[60]*tmpObjS[58] + tmpFx[74]*tmpObjS[71] + tmpFx[88]*tmpObjS[84] + tmpFx[102]*tmpObjS[97] + tmpFx[116]*tmpObjS[110] + tmpFx[130]*tmpObjS[123] + tmpFx[144]*tmpObjS[136] + tmpFx[158]*tmpObjS[149] + tmpFx[172]*tmpObjS[162];
tmpQ2[59] = + tmpFx[4]*tmpObjS[7] + tmpFx[18]*tmpObjS[20] + tmpFx[32]*tmpObjS[33] + tmpFx[46]*tmpObjS[46] + tmpFx[60]*tmpObjS[59] + tmpFx[74]*tmpObjS[72] + tmpFx[88]*tmpObjS[85] + tmpFx[102]*tmpObjS[98] + tmpFx[116]*tmpObjS[111] + tmpFx[130]*tmpObjS[124] + tmpFx[144]*tmpObjS[137] + tmpFx[158]*tmpObjS[150] + tmpFx[172]*tmpObjS[163];
tmpQ2[60] = + tmpFx[4]*tmpObjS[8] + tmpFx[18]*tmpObjS[21] + tmpFx[32]*tmpObjS[34] + tmpFx[46]*tmpObjS[47] + tmpFx[60]*tmpObjS[60] + tmpFx[74]*tmpObjS[73] + tmpFx[88]*tmpObjS[86] + tmpFx[102]*tmpObjS[99] + tmpFx[116]*tmpObjS[112] + tmpFx[130]*tmpObjS[125] + tmpFx[144]*tmpObjS[138] + tmpFx[158]*tmpObjS[151] + tmpFx[172]*tmpObjS[164];
tmpQ2[61] = + tmpFx[4]*tmpObjS[9] + tmpFx[18]*tmpObjS[22] + tmpFx[32]*tmpObjS[35] + tmpFx[46]*tmpObjS[48] + tmpFx[60]*tmpObjS[61] + tmpFx[74]*tmpObjS[74] + tmpFx[88]*tmpObjS[87] + tmpFx[102]*tmpObjS[100] + tmpFx[116]*tmpObjS[113] + tmpFx[130]*tmpObjS[126] + tmpFx[144]*tmpObjS[139] + tmpFx[158]*tmpObjS[152] + tmpFx[172]*tmpObjS[165];
tmpQ2[62] = + tmpFx[4]*tmpObjS[10] + tmpFx[18]*tmpObjS[23] + tmpFx[32]*tmpObjS[36] + tmpFx[46]*tmpObjS[49] + tmpFx[60]*tmpObjS[62] + tmpFx[74]*tmpObjS[75] + tmpFx[88]*tmpObjS[88] + tmpFx[102]*tmpObjS[101] + tmpFx[116]*tmpObjS[114] + tmpFx[130]*tmpObjS[127] + tmpFx[144]*tmpObjS[140] + tmpFx[158]*tmpObjS[153] + tmpFx[172]*tmpObjS[166];
tmpQ2[63] = + tmpFx[4]*tmpObjS[11] + tmpFx[18]*tmpObjS[24] + tmpFx[32]*tmpObjS[37] + tmpFx[46]*tmpObjS[50] + tmpFx[60]*tmpObjS[63] + tmpFx[74]*tmpObjS[76] + tmpFx[88]*tmpObjS[89] + tmpFx[102]*tmpObjS[102] + tmpFx[116]*tmpObjS[115] + tmpFx[130]*tmpObjS[128] + tmpFx[144]*tmpObjS[141] + tmpFx[158]*tmpObjS[154] + tmpFx[172]*tmpObjS[167];
tmpQ2[64] = + tmpFx[4]*tmpObjS[12] + tmpFx[18]*tmpObjS[25] + tmpFx[32]*tmpObjS[38] + tmpFx[46]*tmpObjS[51] + tmpFx[60]*tmpObjS[64] + tmpFx[74]*tmpObjS[77] + tmpFx[88]*tmpObjS[90] + tmpFx[102]*tmpObjS[103] + tmpFx[116]*tmpObjS[116] + tmpFx[130]*tmpObjS[129] + tmpFx[144]*tmpObjS[142] + tmpFx[158]*tmpObjS[155] + tmpFx[172]*tmpObjS[168];
tmpQ2[65] = + tmpFx[5]*tmpObjS[0] + tmpFx[19]*tmpObjS[13] + tmpFx[33]*tmpObjS[26] + tmpFx[47]*tmpObjS[39] + tmpFx[61]*tmpObjS[52] + tmpFx[75]*tmpObjS[65] + tmpFx[89]*tmpObjS[78] + tmpFx[103]*tmpObjS[91] + tmpFx[117]*tmpObjS[104] + tmpFx[131]*tmpObjS[117] + tmpFx[145]*tmpObjS[130] + tmpFx[159]*tmpObjS[143] + tmpFx[173]*tmpObjS[156];
tmpQ2[66] = + tmpFx[5]*tmpObjS[1] + tmpFx[19]*tmpObjS[14] + tmpFx[33]*tmpObjS[27] + tmpFx[47]*tmpObjS[40] + tmpFx[61]*tmpObjS[53] + tmpFx[75]*tmpObjS[66] + tmpFx[89]*tmpObjS[79] + tmpFx[103]*tmpObjS[92] + tmpFx[117]*tmpObjS[105] + tmpFx[131]*tmpObjS[118] + tmpFx[145]*tmpObjS[131] + tmpFx[159]*tmpObjS[144] + tmpFx[173]*tmpObjS[157];
tmpQ2[67] = + tmpFx[5]*tmpObjS[2] + tmpFx[19]*tmpObjS[15] + tmpFx[33]*tmpObjS[28] + tmpFx[47]*tmpObjS[41] + tmpFx[61]*tmpObjS[54] + tmpFx[75]*tmpObjS[67] + tmpFx[89]*tmpObjS[80] + tmpFx[103]*tmpObjS[93] + tmpFx[117]*tmpObjS[106] + tmpFx[131]*tmpObjS[119] + tmpFx[145]*tmpObjS[132] + tmpFx[159]*tmpObjS[145] + tmpFx[173]*tmpObjS[158];
tmpQ2[68] = + tmpFx[5]*tmpObjS[3] + tmpFx[19]*tmpObjS[16] + tmpFx[33]*tmpObjS[29] + tmpFx[47]*tmpObjS[42] + tmpFx[61]*tmpObjS[55] + tmpFx[75]*tmpObjS[68] + tmpFx[89]*tmpObjS[81] + tmpFx[103]*tmpObjS[94] + tmpFx[117]*tmpObjS[107] + tmpFx[131]*tmpObjS[120] + tmpFx[145]*tmpObjS[133] + tmpFx[159]*tmpObjS[146] + tmpFx[173]*tmpObjS[159];
tmpQ2[69] = + tmpFx[5]*tmpObjS[4] + tmpFx[19]*tmpObjS[17] + tmpFx[33]*tmpObjS[30] + tmpFx[47]*tmpObjS[43] + tmpFx[61]*tmpObjS[56] + tmpFx[75]*tmpObjS[69] + tmpFx[89]*tmpObjS[82] + tmpFx[103]*tmpObjS[95] + tmpFx[117]*tmpObjS[108] + tmpFx[131]*tmpObjS[121] + tmpFx[145]*tmpObjS[134] + tmpFx[159]*tmpObjS[147] + tmpFx[173]*tmpObjS[160];
tmpQ2[70] = + tmpFx[5]*tmpObjS[5] + tmpFx[19]*tmpObjS[18] + tmpFx[33]*tmpObjS[31] + tmpFx[47]*tmpObjS[44] + tmpFx[61]*tmpObjS[57] + tmpFx[75]*tmpObjS[70] + tmpFx[89]*tmpObjS[83] + tmpFx[103]*tmpObjS[96] + tmpFx[117]*tmpObjS[109] + tmpFx[131]*tmpObjS[122] + tmpFx[145]*tmpObjS[135] + tmpFx[159]*tmpObjS[148] + tmpFx[173]*tmpObjS[161];
tmpQ2[71] = + tmpFx[5]*tmpObjS[6] + tmpFx[19]*tmpObjS[19] + tmpFx[33]*tmpObjS[32] + tmpFx[47]*tmpObjS[45] + tmpFx[61]*tmpObjS[58] + tmpFx[75]*tmpObjS[71] + tmpFx[89]*tmpObjS[84] + tmpFx[103]*tmpObjS[97] + tmpFx[117]*tmpObjS[110] + tmpFx[131]*tmpObjS[123] + tmpFx[145]*tmpObjS[136] + tmpFx[159]*tmpObjS[149] + tmpFx[173]*tmpObjS[162];
tmpQ2[72] = + tmpFx[5]*tmpObjS[7] + tmpFx[19]*tmpObjS[20] + tmpFx[33]*tmpObjS[33] + tmpFx[47]*tmpObjS[46] + tmpFx[61]*tmpObjS[59] + tmpFx[75]*tmpObjS[72] + tmpFx[89]*tmpObjS[85] + tmpFx[103]*tmpObjS[98] + tmpFx[117]*tmpObjS[111] + tmpFx[131]*tmpObjS[124] + tmpFx[145]*tmpObjS[137] + tmpFx[159]*tmpObjS[150] + tmpFx[173]*tmpObjS[163];
tmpQ2[73] = + tmpFx[5]*tmpObjS[8] + tmpFx[19]*tmpObjS[21] + tmpFx[33]*tmpObjS[34] + tmpFx[47]*tmpObjS[47] + tmpFx[61]*tmpObjS[60] + tmpFx[75]*tmpObjS[73] + tmpFx[89]*tmpObjS[86] + tmpFx[103]*tmpObjS[99] + tmpFx[117]*tmpObjS[112] + tmpFx[131]*tmpObjS[125] + tmpFx[145]*tmpObjS[138] + tmpFx[159]*tmpObjS[151] + tmpFx[173]*tmpObjS[164];
tmpQ2[74] = + tmpFx[5]*tmpObjS[9] + tmpFx[19]*tmpObjS[22] + tmpFx[33]*tmpObjS[35] + tmpFx[47]*tmpObjS[48] + tmpFx[61]*tmpObjS[61] + tmpFx[75]*tmpObjS[74] + tmpFx[89]*tmpObjS[87] + tmpFx[103]*tmpObjS[100] + tmpFx[117]*tmpObjS[113] + tmpFx[131]*tmpObjS[126] + tmpFx[145]*tmpObjS[139] + tmpFx[159]*tmpObjS[152] + tmpFx[173]*tmpObjS[165];
tmpQ2[75] = + tmpFx[5]*tmpObjS[10] + tmpFx[19]*tmpObjS[23] + tmpFx[33]*tmpObjS[36] + tmpFx[47]*tmpObjS[49] + tmpFx[61]*tmpObjS[62] + tmpFx[75]*tmpObjS[75] + tmpFx[89]*tmpObjS[88] + tmpFx[103]*tmpObjS[101] + tmpFx[117]*tmpObjS[114] + tmpFx[131]*tmpObjS[127] + tmpFx[145]*tmpObjS[140] + tmpFx[159]*tmpObjS[153] + tmpFx[173]*tmpObjS[166];
tmpQ2[76] = + tmpFx[5]*tmpObjS[11] + tmpFx[19]*tmpObjS[24] + tmpFx[33]*tmpObjS[37] + tmpFx[47]*tmpObjS[50] + tmpFx[61]*tmpObjS[63] + tmpFx[75]*tmpObjS[76] + tmpFx[89]*tmpObjS[89] + tmpFx[103]*tmpObjS[102] + tmpFx[117]*tmpObjS[115] + tmpFx[131]*tmpObjS[128] + tmpFx[145]*tmpObjS[141] + tmpFx[159]*tmpObjS[154] + tmpFx[173]*tmpObjS[167];
tmpQ2[77] = + tmpFx[5]*tmpObjS[12] + tmpFx[19]*tmpObjS[25] + tmpFx[33]*tmpObjS[38] + tmpFx[47]*tmpObjS[51] + tmpFx[61]*tmpObjS[64] + tmpFx[75]*tmpObjS[77] + tmpFx[89]*tmpObjS[90] + tmpFx[103]*tmpObjS[103] + tmpFx[117]*tmpObjS[116] + tmpFx[131]*tmpObjS[129] + tmpFx[145]*tmpObjS[142] + tmpFx[159]*tmpObjS[155] + tmpFx[173]*tmpObjS[168];
tmpQ2[78] = + tmpFx[6]*tmpObjS[0] + tmpFx[20]*tmpObjS[13] + tmpFx[34]*tmpObjS[26] + tmpFx[48]*tmpObjS[39] + tmpFx[62]*tmpObjS[52] + tmpFx[76]*tmpObjS[65] + tmpFx[90]*tmpObjS[78] + tmpFx[104]*tmpObjS[91] + tmpFx[118]*tmpObjS[104] + tmpFx[132]*tmpObjS[117] + tmpFx[146]*tmpObjS[130] + tmpFx[160]*tmpObjS[143] + tmpFx[174]*tmpObjS[156];
tmpQ2[79] = + tmpFx[6]*tmpObjS[1] + tmpFx[20]*tmpObjS[14] + tmpFx[34]*tmpObjS[27] + tmpFx[48]*tmpObjS[40] + tmpFx[62]*tmpObjS[53] + tmpFx[76]*tmpObjS[66] + tmpFx[90]*tmpObjS[79] + tmpFx[104]*tmpObjS[92] + tmpFx[118]*tmpObjS[105] + tmpFx[132]*tmpObjS[118] + tmpFx[146]*tmpObjS[131] + tmpFx[160]*tmpObjS[144] + tmpFx[174]*tmpObjS[157];
tmpQ2[80] = + tmpFx[6]*tmpObjS[2] + tmpFx[20]*tmpObjS[15] + tmpFx[34]*tmpObjS[28] + tmpFx[48]*tmpObjS[41] + tmpFx[62]*tmpObjS[54] + tmpFx[76]*tmpObjS[67] + tmpFx[90]*tmpObjS[80] + tmpFx[104]*tmpObjS[93] + tmpFx[118]*tmpObjS[106] + tmpFx[132]*tmpObjS[119] + tmpFx[146]*tmpObjS[132] + tmpFx[160]*tmpObjS[145] + tmpFx[174]*tmpObjS[158];
tmpQ2[81] = + tmpFx[6]*tmpObjS[3] + tmpFx[20]*tmpObjS[16] + tmpFx[34]*tmpObjS[29] + tmpFx[48]*tmpObjS[42] + tmpFx[62]*tmpObjS[55] + tmpFx[76]*tmpObjS[68] + tmpFx[90]*tmpObjS[81] + tmpFx[104]*tmpObjS[94] + tmpFx[118]*tmpObjS[107] + tmpFx[132]*tmpObjS[120] + tmpFx[146]*tmpObjS[133] + tmpFx[160]*tmpObjS[146] + tmpFx[174]*tmpObjS[159];
tmpQ2[82] = + tmpFx[6]*tmpObjS[4] + tmpFx[20]*tmpObjS[17] + tmpFx[34]*tmpObjS[30] + tmpFx[48]*tmpObjS[43] + tmpFx[62]*tmpObjS[56] + tmpFx[76]*tmpObjS[69] + tmpFx[90]*tmpObjS[82] + tmpFx[104]*tmpObjS[95] + tmpFx[118]*tmpObjS[108] + tmpFx[132]*tmpObjS[121] + tmpFx[146]*tmpObjS[134] + tmpFx[160]*tmpObjS[147] + tmpFx[174]*tmpObjS[160];
tmpQ2[83] = + tmpFx[6]*tmpObjS[5] + tmpFx[20]*tmpObjS[18] + tmpFx[34]*tmpObjS[31] + tmpFx[48]*tmpObjS[44] + tmpFx[62]*tmpObjS[57] + tmpFx[76]*tmpObjS[70] + tmpFx[90]*tmpObjS[83] + tmpFx[104]*tmpObjS[96] + tmpFx[118]*tmpObjS[109] + tmpFx[132]*tmpObjS[122] + tmpFx[146]*tmpObjS[135] + tmpFx[160]*tmpObjS[148] + tmpFx[174]*tmpObjS[161];
tmpQ2[84] = + tmpFx[6]*tmpObjS[6] + tmpFx[20]*tmpObjS[19] + tmpFx[34]*tmpObjS[32] + tmpFx[48]*tmpObjS[45] + tmpFx[62]*tmpObjS[58] + tmpFx[76]*tmpObjS[71] + tmpFx[90]*tmpObjS[84] + tmpFx[104]*tmpObjS[97] + tmpFx[118]*tmpObjS[110] + tmpFx[132]*tmpObjS[123] + tmpFx[146]*tmpObjS[136] + tmpFx[160]*tmpObjS[149] + tmpFx[174]*tmpObjS[162];
tmpQ2[85] = + tmpFx[6]*tmpObjS[7] + tmpFx[20]*tmpObjS[20] + tmpFx[34]*tmpObjS[33] + tmpFx[48]*tmpObjS[46] + tmpFx[62]*tmpObjS[59] + tmpFx[76]*tmpObjS[72] + tmpFx[90]*tmpObjS[85] + tmpFx[104]*tmpObjS[98] + tmpFx[118]*tmpObjS[111] + tmpFx[132]*tmpObjS[124] + tmpFx[146]*tmpObjS[137] + tmpFx[160]*tmpObjS[150] + tmpFx[174]*tmpObjS[163];
tmpQ2[86] = + tmpFx[6]*tmpObjS[8] + tmpFx[20]*tmpObjS[21] + tmpFx[34]*tmpObjS[34] + tmpFx[48]*tmpObjS[47] + tmpFx[62]*tmpObjS[60] + tmpFx[76]*tmpObjS[73] + tmpFx[90]*tmpObjS[86] + tmpFx[104]*tmpObjS[99] + tmpFx[118]*tmpObjS[112] + tmpFx[132]*tmpObjS[125] + tmpFx[146]*tmpObjS[138] + tmpFx[160]*tmpObjS[151] + tmpFx[174]*tmpObjS[164];
tmpQ2[87] = + tmpFx[6]*tmpObjS[9] + tmpFx[20]*tmpObjS[22] + tmpFx[34]*tmpObjS[35] + tmpFx[48]*tmpObjS[48] + tmpFx[62]*tmpObjS[61] + tmpFx[76]*tmpObjS[74] + tmpFx[90]*tmpObjS[87] + tmpFx[104]*tmpObjS[100] + tmpFx[118]*tmpObjS[113] + tmpFx[132]*tmpObjS[126] + tmpFx[146]*tmpObjS[139] + tmpFx[160]*tmpObjS[152] + tmpFx[174]*tmpObjS[165];
tmpQ2[88] = + tmpFx[6]*tmpObjS[10] + tmpFx[20]*tmpObjS[23] + tmpFx[34]*tmpObjS[36] + tmpFx[48]*tmpObjS[49] + tmpFx[62]*tmpObjS[62] + tmpFx[76]*tmpObjS[75] + tmpFx[90]*tmpObjS[88] + tmpFx[104]*tmpObjS[101] + tmpFx[118]*tmpObjS[114] + tmpFx[132]*tmpObjS[127] + tmpFx[146]*tmpObjS[140] + tmpFx[160]*tmpObjS[153] + tmpFx[174]*tmpObjS[166];
tmpQ2[89] = + tmpFx[6]*tmpObjS[11] + tmpFx[20]*tmpObjS[24] + tmpFx[34]*tmpObjS[37] + tmpFx[48]*tmpObjS[50] + tmpFx[62]*tmpObjS[63] + tmpFx[76]*tmpObjS[76] + tmpFx[90]*tmpObjS[89] + tmpFx[104]*tmpObjS[102] + tmpFx[118]*tmpObjS[115] + tmpFx[132]*tmpObjS[128] + tmpFx[146]*tmpObjS[141] + tmpFx[160]*tmpObjS[154] + tmpFx[174]*tmpObjS[167];
tmpQ2[90] = + tmpFx[6]*tmpObjS[12] + tmpFx[20]*tmpObjS[25] + tmpFx[34]*tmpObjS[38] + tmpFx[48]*tmpObjS[51] + tmpFx[62]*tmpObjS[64] + tmpFx[76]*tmpObjS[77] + tmpFx[90]*tmpObjS[90] + tmpFx[104]*tmpObjS[103] + tmpFx[118]*tmpObjS[116] + tmpFx[132]*tmpObjS[129] + tmpFx[146]*tmpObjS[142] + tmpFx[160]*tmpObjS[155] + tmpFx[174]*tmpObjS[168];
tmpQ2[91] = + tmpFx[7]*tmpObjS[0] + tmpFx[21]*tmpObjS[13] + tmpFx[35]*tmpObjS[26] + tmpFx[49]*tmpObjS[39] + tmpFx[63]*tmpObjS[52] + tmpFx[77]*tmpObjS[65] + tmpFx[91]*tmpObjS[78] + tmpFx[105]*tmpObjS[91] + tmpFx[119]*tmpObjS[104] + tmpFx[133]*tmpObjS[117] + tmpFx[147]*tmpObjS[130] + tmpFx[161]*tmpObjS[143] + tmpFx[175]*tmpObjS[156];
tmpQ2[92] = + tmpFx[7]*tmpObjS[1] + tmpFx[21]*tmpObjS[14] + tmpFx[35]*tmpObjS[27] + tmpFx[49]*tmpObjS[40] + tmpFx[63]*tmpObjS[53] + tmpFx[77]*tmpObjS[66] + tmpFx[91]*tmpObjS[79] + tmpFx[105]*tmpObjS[92] + tmpFx[119]*tmpObjS[105] + tmpFx[133]*tmpObjS[118] + tmpFx[147]*tmpObjS[131] + tmpFx[161]*tmpObjS[144] + tmpFx[175]*tmpObjS[157];
tmpQ2[93] = + tmpFx[7]*tmpObjS[2] + tmpFx[21]*tmpObjS[15] + tmpFx[35]*tmpObjS[28] + tmpFx[49]*tmpObjS[41] + tmpFx[63]*tmpObjS[54] + tmpFx[77]*tmpObjS[67] + tmpFx[91]*tmpObjS[80] + tmpFx[105]*tmpObjS[93] + tmpFx[119]*tmpObjS[106] + tmpFx[133]*tmpObjS[119] + tmpFx[147]*tmpObjS[132] + tmpFx[161]*tmpObjS[145] + tmpFx[175]*tmpObjS[158];
tmpQ2[94] = + tmpFx[7]*tmpObjS[3] + tmpFx[21]*tmpObjS[16] + tmpFx[35]*tmpObjS[29] + tmpFx[49]*tmpObjS[42] + tmpFx[63]*tmpObjS[55] + tmpFx[77]*tmpObjS[68] + tmpFx[91]*tmpObjS[81] + tmpFx[105]*tmpObjS[94] + tmpFx[119]*tmpObjS[107] + tmpFx[133]*tmpObjS[120] + tmpFx[147]*tmpObjS[133] + tmpFx[161]*tmpObjS[146] + tmpFx[175]*tmpObjS[159];
tmpQ2[95] = + tmpFx[7]*tmpObjS[4] + tmpFx[21]*tmpObjS[17] + tmpFx[35]*tmpObjS[30] + tmpFx[49]*tmpObjS[43] + tmpFx[63]*tmpObjS[56] + tmpFx[77]*tmpObjS[69] + tmpFx[91]*tmpObjS[82] + tmpFx[105]*tmpObjS[95] + tmpFx[119]*tmpObjS[108] + tmpFx[133]*tmpObjS[121] + tmpFx[147]*tmpObjS[134] + tmpFx[161]*tmpObjS[147] + tmpFx[175]*tmpObjS[160];
tmpQ2[96] = + tmpFx[7]*tmpObjS[5] + tmpFx[21]*tmpObjS[18] + tmpFx[35]*tmpObjS[31] + tmpFx[49]*tmpObjS[44] + tmpFx[63]*tmpObjS[57] + tmpFx[77]*tmpObjS[70] + tmpFx[91]*tmpObjS[83] + tmpFx[105]*tmpObjS[96] + tmpFx[119]*tmpObjS[109] + tmpFx[133]*tmpObjS[122] + tmpFx[147]*tmpObjS[135] + tmpFx[161]*tmpObjS[148] + tmpFx[175]*tmpObjS[161];
tmpQ2[97] = + tmpFx[7]*tmpObjS[6] + tmpFx[21]*tmpObjS[19] + tmpFx[35]*tmpObjS[32] + tmpFx[49]*tmpObjS[45] + tmpFx[63]*tmpObjS[58] + tmpFx[77]*tmpObjS[71] + tmpFx[91]*tmpObjS[84] + tmpFx[105]*tmpObjS[97] + tmpFx[119]*tmpObjS[110] + tmpFx[133]*tmpObjS[123] + tmpFx[147]*tmpObjS[136] + tmpFx[161]*tmpObjS[149] + tmpFx[175]*tmpObjS[162];
tmpQ2[98] = + tmpFx[7]*tmpObjS[7] + tmpFx[21]*tmpObjS[20] + tmpFx[35]*tmpObjS[33] + tmpFx[49]*tmpObjS[46] + tmpFx[63]*tmpObjS[59] + tmpFx[77]*tmpObjS[72] + tmpFx[91]*tmpObjS[85] + tmpFx[105]*tmpObjS[98] + tmpFx[119]*tmpObjS[111] + tmpFx[133]*tmpObjS[124] + tmpFx[147]*tmpObjS[137] + tmpFx[161]*tmpObjS[150] + tmpFx[175]*tmpObjS[163];
tmpQ2[99] = + tmpFx[7]*tmpObjS[8] + tmpFx[21]*tmpObjS[21] + tmpFx[35]*tmpObjS[34] + tmpFx[49]*tmpObjS[47] + tmpFx[63]*tmpObjS[60] + tmpFx[77]*tmpObjS[73] + tmpFx[91]*tmpObjS[86] + tmpFx[105]*tmpObjS[99] + tmpFx[119]*tmpObjS[112] + tmpFx[133]*tmpObjS[125] + tmpFx[147]*tmpObjS[138] + tmpFx[161]*tmpObjS[151] + tmpFx[175]*tmpObjS[164];
tmpQ2[100] = + tmpFx[7]*tmpObjS[9] + tmpFx[21]*tmpObjS[22] + tmpFx[35]*tmpObjS[35] + tmpFx[49]*tmpObjS[48] + tmpFx[63]*tmpObjS[61] + tmpFx[77]*tmpObjS[74] + tmpFx[91]*tmpObjS[87] + tmpFx[105]*tmpObjS[100] + tmpFx[119]*tmpObjS[113] + tmpFx[133]*tmpObjS[126] + tmpFx[147]*tmpObjS[139] + tmpFx[161]*tmpObjS[152] + tmpFx[175]*tmpObjS[165];
tmpQ2[101] = + tmpFx[7]*tmpObjS[10] + tmpFx[21]*tmpObjS[23] + tmpFx[35]*tmpObjS[36] + tmpFx[49]*tmpObjS[49] + tmpFx[63]*tmpObjS[62] + tmpFx[77]*tmpObjS[75] + tmpFx[91]*tmpObjS[88] + tmpFx[105]*tmpObjS[101] + tmpFx[119]*tmpObjS[114] + tmpFx[133]*tmpObjS[127] + tmpFx[147]*tmpObjS[140] + tmpFx[161]*tmpObjS[153] + tmpFx[175]*tmpObjS[166];
tmpQ2[102] = + tmpFx[7]*tmpObjS[11] + tmpFx[21]*tmpObjS[24] + tmpFx[35]*tmpObjS[37] + tmpFx[49]*tmpObjS[50] + tmpFx[63]*tmpObjS[63] + tmpFx[77]*tmpObjS[76] + tmpFx[91]*tmpObjS[89] + tmpFx[105]*tmpObjS[102] + tmpFx[119]*tmpObjS[115] + tmpFx[133]*tmpObjS[128] + tmpFx[147]*tmpObjS[141] + tmpFx[161]*tmpObjS[154] + tmpFx[175]*tmpObjS[167];
tmpQ2[103] = + tmpFx[7]*tmpObjS[12] + tmpFx[21]*tmpObjS[25] + tmpFx[35]*tmpObjS[38] + tmpFx[49]*tmpObjS[51] + tmpFx[63]*tmpObjS[64] + tmpFx[77]*tmpObjS[77] + tmpFx[91]*tmpObjS[90] + tmpFx[105]*tmpObjS[103] + tmpFx[119]*tmpObjS[116] + tmpFx[133]*tmpObjS[129] + tmpFx[147]*tmpObjS[142] + tmpFx[161]*tmpObjS[155] + tmpFx[175]*tmpObjS[168];
tmpQ2[104] = + tmpFx[8]*tmpObjS[0] + tmpFx[22]*tmpObjS[13] + tmpFx[36]*tmpObjS[26] + tmpFx[50]*tmpObjS[39] + tmpFx[64]*tmpObjS[52] + tmpFx[78]*tmpObjS[65] + tmpFx[92]*tmpObjS[78] + tmpFx[106]*tmpObjS[91] + tmpFx[120]*tmpObjS[104] + tmpFx[134]*tmpObjS[117] + tmpFx[148]*tmpObjS[130] + tmpFx[162]*tmpObjS[143] + tmpFx[176]*tmpObjS[156];
tmpQ2[105] = + tmpFx[8]*tmpObjS[1] + tmpFx[22]*tmpObjS[14] + tmpFx[36]*tmpObjS[27] + tmpFx[50]*tmpObjS[40] + tmpFx[64]*tmpObjS[53] + tmpFx[78]*tmpObjS[66] + tmpFx[92]*tmpObjS[79] + tmpFx[106]*tmpObjS[92] + tmpFx[120]*tmpObjS[105] + tmpFx[134]*tmpObjS[118] + tmpFx[148]*tmpObjS[131] + tmpFx[162]*tmpObjS[144] + tmpFx[176]*tmpObjS[157];
tmpQ2[106] = + tmpFx[8]*tmpObjS[2] + tmpFx[22]*tmpObjS[15] + tmpFx[36]*tmpObjS[28] + tmpFx[50]*tmpObjS[41] + tmpFx[64]*tmpObjS[54] + tmpFx[78]*tmpObjS[67] + tmpFx[92]*tmpObjS[80] + tmpFx[106]*tmpObjS[93] + tmpFx[120]*tmpObjS[106] + tmpFx[134]*tmpObjS[119] + tmpFx[148]*tmpObjS[132] + tmpFx[162]*tmpObjS[145] + tmpFx[176]*tmpObjS[158];
tmpQ2[107] = + tmpFx[8]*tmpObjS[3] + tmpFx[22]*tmpObjS[16] + tmpFx[36]*tmpObjS[29] + tmpFx[50]*tmpObjS[42] + tmpFx[64]*tmpObjS[55] + tmpFx[78]*tmpObjS[68] + tmpFx[92]*tmpObjS[81] + tmpFx[106]*tmpObjS[94] + tmpFx[120]*tmpObjS[107] + tmpFx[134]*tmpObjS[120] + tmpFx[148]*tmpObjS[133] + tmpFx[162]*tmpObjS[146] + tmpFx[176]*tmpObjS[159];
tmpQ2[108] = + tmpFx[8]*tmpObjS[4] + tmpFx[22]*tmpObjS[17] + tmpFx[36]*tmpObjS[30] + tmpFx[50]*tmpObjS[43] + tmpFx[64]*tmpObjS[56] + tmpFx[78]*tmpObjS[69] + tmpFx[92]*tmpObjS[82] + tmpFx[106]*tmpObjS[95] + tmpFx[120]*tmpObjS[108] + tmpFx[134]*tmpObjS[121] + tmpFx[148]*tmpObjS[134] + tmpFx[162]*tmpObjS[147] + tmpFx[176]*tmpObjS[160];
tmpQ2[109] = + tmpFx[8]*tmpObjS[5] + tmpFx[22]*tmpObjS[18] + tmpFx[36]*tmpObjS[31] + tmpFx[50]*tmpObjS[44] + tmpFx[64]*tmpObjS[57] + tmpFx[78]*tmpObjS[70] + tmpFx[92]*tmpObjS[83] + tmpFx[106]*tmpObjS[96] + tmpFx[120]*tmpObjS[109] + tmpFx[134]*tmpObjS[122] + tmpFx[148]*tmpObjS[135] + tmpFx[162]*tmpObjS[148] + tmpFx[176]*tmpObjS[161];
tmpQ2[110] = + tmpFx[8]*tmpObjS[6] + tmpFx[22]*tmpObjS[19] + tmpFx[36]*tmpObjS[32] + tmpFx[50]*tmpObjS[45] + tmpFx[64]*tmpObjS[58] + tmpFx[78]*tmpObjS[71] + tmpFx[92]*tmpObjS[84] + tmpFx[106]*tmpObjS[97] + tmpFx[120]*tmpObjS[110] + tmpFx[134]*tmpObjS[123] + tmpFx[148]*tmpObjS[136] + tmpFx[162]*tmpObjS[149] + tmpFx[176]*tmpObjS[162];
tmpQ2[111] = + tmpFx[8]*tmpObjS[7] + tmpFx[22]*tmpObjS[20] + tmpFx[36]*tmpObjS[33] + tmpFx[50]*tmpObjS[46] + tmpFx[64]*tmpObjS[59] + tmpFx[78]*tmpObjS[72] + tmpFx[92]*tmpObjS[85] + tmpFx[106]*tmpObjS[98] + tmpFx[120]*tmpObjS[111] + tmpFx[134]*tmpObjS[124] + tmpFx[148]*tmpObjS[137] + tmpFx[162]*tmpObjS[150] + tmpFx[176]*tmpObjS[163];
tmpQ2[112] = + tmpFx[8]*tmpObjS[8] + tmpFx[22]*tmpObjS[21] + tmpFx[36]*tmpObjS[34] + tmpFx[50]*tmpObjS[47] + tmpFx[64]*tmpObjS[60] + tmpFx[78]*tmpObjS[73] + tmpFx[92]*tmpObjS[86] + tmpFx[106]*tmpObjS[99] + tmpFx[120]*tmpObjS[112] + tmpFx[134]*tmpObjS[125] + tmpFx[148]*tmpObjS[138] + tmpFx[162]*tmpObjS[151] + tmpFx[176]*tmpObjS[164];
tmpQ2[113] = + tmpFx[8]*tmpObjS[9] + tmpFx[22]*tmpObjS[22] + tmpFx[36]*tmpObjS[35] + tmpFx[50]*tmpObjS[48] + tmpFx[64]*tmpObjS[61] + tmpFx[78]*tmpObjS[74] + tmpFx[92]*tmpObjS[87] + tmpFx[106]*tmpObjS[100] + tmpFx[120]*tmpObjS[113] + tmpFx[134]*tmpObjS[126] + tmpFx[148]*tmpObjS[139] + tmpFx[162]*tmpObjS[152] + tmpFx[176]*tmpObjS[165];
tmpQ2[114] = + tmpFx[8]*tmpObjS[10] + tmpFx[22]*tmpObjS[23] + tmpFx[36]*tmpObjS[36] + tmpFx[50]*tmpObjS[49] + tmpFx[64]*tmpObjS[62] + tmpFx[78]*tmpObjS[75] + tmpFx[92]*tmpObjS[88] + tmpFx[106]*tmpObjS[101] + tmpFx[120]*tmpObjS[114] + tmpFx[134]*tmpObjS[127] + tmpFx[148]*tmpObjS[140] + tmpFx[162]*tmpObjS[153] + tmpFx[176]*tmpObjS[166];
tmpQ2[115] = + tmpFx[8]*tmpObjS[11] + tmpFx[22]*tmpObjS[24] + tmpFx[36]*tmpObjS[37] + tmpFx[50]*tmpObjS[50] + tmpFx[64]*tmpObjS[63] + tmpFx[78]*tmpObjS[76] + tmpFx[92]*tmpObjS[89] + tmpFx[106]*tmpObjS[102] + tmpFx[120]*tmpObjS[115] + tmpFx[134]*tmpObjS[128] + tmpFx[148]*tmpObjS[141] + tmpFx[162]*tmpObjS[154] + tmpFx[176]*tmpObjS[167];
tmpQ2[116] = + tmpFx[8]*tmpObjS[12] + tmpFx[22]*tmpObjS[25] + tmpFx[36]*tmpObjS[38] + tmpFx[50]*tmpObjS[51] + tmpFx[64]*tmpObjS[64] + tmpFx[78]*tmpObjS[77] + tmpFx[92]*tmpObjS[90] + tmpFx[106]*tmpObjS[103] + tmpFx[120]*tmpObjS[116] + tmpFx[134]*tmpObjS[129] + tmpFx[148]*tmpObjS[142] + tmpFx[162]*tmpObjS[155] + tmpFx[176]*tmpObjS[168];
tmpQ2[117] = + tmpFx[9]*tmpObjS[0] + tmpFx[23]*tmpObjS[13] + tmpFx[37]*tmpObjS[26] + tmpFx[51]*tmpObjS[39] + tmpFx[65]*tmpObjS[52] + tmpFx[79]*tmpObjS[65] + tmpFx[93]*tmpObjS[78] + tmpFx[107]*tmpObjS[91] + tmpFx[121]*tmpObjS[104] + tmpFx[135]*tmpObjS[117] + tmpFx[149]*tmpObjS[130] + tmpFx[163]*tmpObjS[143] + tmpFx[177]*tmpObjS[156];
tmpQ2[118] = + tmpFx[9]*tmpObjS[1] + tmpFx[23]*tmpObjS[14] + tmpFx[37]*tmpObjS[27] + tmpFx[51]*tmpObjS[40] + tmpFx[65]*tmpObjS[53] + tmpFx[79]*tmpObjS[66] + tmpFx[93]*tmpObjS[79] + tmpFx[107]*tmpObjS[92] + tmpFx[121]*tmpObjS[105] + tmpFx[135]*tmpObjS[118] + tmpFx[149]*tmpObjS[131] + tmpFx[163]*tmpObjS[144] + tmpFx[177]*tmpObjS[157];
tmpQ2[119] = + tmpFx[9]*tmpObjS[2] + tmpFx[23]*tmpObjS[15] + tmpFx[37]*tmpObjS[28] + tmpFx[51]*tmpObjS[41] + tmpFx[65]*tmpObjS[54] + tmpFx[79]*tmpObjS[67] + tmpFx[93]*tmpObjS[80] + tmpFx[107]*tmpObjS[93] + tmpFx[121]*tmpObjS[106] + tmpFx[135]*tmpObjS[119] + tmpFx[149]*tmpObjS[132] + tmpFx[163]*tmpObjS[145] + tmpFx[177]*tmpObjS[158];
tmpQ2[120] = + tmpFx[9]*tmpObjS[3] + tmpFx[23]*tmpObjS[16] + tmpFx[37]*tmpObjS[29] + tmpFx[51]*tmpObjS[42] + tmpFx[65]*tmpObjS[55] + tmpFx[79]*tmpObjS[68] + tmpFx[93]*tmpObjS[81] + tmpFx[107]*tmpObjS[94] + tmpFx[121]*tmpObjS[107] + tmpFx[135]*tmpObjS[120] + tmpFx[149]*tmpObjS[133] + tmpFx[163]*tmpObjS[146] + tmpFx[177]*tmpObjS[159];
tmpQ2[121] = + tmpFx[9]*tmpObjS[4] + tmpFx[23]*tmpObjS[17] + tmpFx[37]*tmpObjS[30] + tmpFx[51]*tmpObjS[43] + tmpFx[65]*tmpObjS[56] + tmpFx[79]*tmpObjS[69] + tmpFx[93]*tmpObjS[82] + tmpFx[107]*tmpObjS[95] + tmpFx[121]*tmpObjS[108] + tmpFx[135]*tmpObjS[121] + tmpFx[149]*tmpObjS[134] + tmpFx[163]*tmpObjS[147] + tmpFx[177]*tmpObjS[160];
tmpQ2[122] = + tmpFx[9]*tmpObjS[5] + tmpFx[23]*tmpObjS[18] + tmpFx[37]*tmpObjS[31] + tmpFx[51]*tmpObjS[44] + tmpFx[65]*tmpObjS[57] + tmpFx[79]*tmpObjS[70] + tmpFx[93]*tmpObjS[83] + tmpFx[107]*tmpObjS[96] + tmpFx[121]*tmpObjS[109] + tmpFx[135]*tmpObjS[122] + tmpFx[149]*tmpObjS[135] + tmpFx[163]*tmpObjS[148] + tmpFx[177]*tmpObjS[161];
tmpQ2[123] = + tmpFx[9]*tmpObjS[6] + tmpFx[23]*tmpObjS[19] + tmpFx[37]*tmpObjS[32] + tmpFx[51]*tmpObjS[45] + tmpFx[65]*tmpObjS[58] + tmpFx[79]*tmpObjS[71] + tmpFx[93]*tmpObjS[84] + tmpFx[107]*tmpObjS[97] + tmpFx[121]*tmpObjS[110] + tmpFx[135]*tmpObjS[123] + tmpFx[149]*tmpObjS[136] + tmpFx[163]*tmpObjS[149] + tmpFx[177]*tmpObjS[162];
tmpQ2[124] = + tmpFx[9]*tmpObjS[7] + tmpFx[23]*tmpObjS[20] + tmpFx[37]*tmpObjS[33] + tmpFx[51]*tmpObjS[46] + tmpFx[65]*tmpObjS[59] + tmpFx[79]*tmpObjS[72] + tmpFx[93]*tmpObjS[85] + tmpFx[107]*tmpObjS[98] + tmpFx[121]*tmpObjS[111] + tmpFx[135]*tmpObjS[124] + tmpFx[149]*tmpObjS[137] + tmpFx[163]*tmpObjS[150] + tmpFx[177]*tmpObjS[163];
tmpQ2[125] = + tmpFx[9]*tmpObjS[8] + tmpFx[23]*tmpObjS[21] + tmpFx[37]*tmpObjS[34] + tmpFx[51]*tmpObjS[47] + tmpFx[65]*tmpObjS[60] + tmpFx[79]*tmpObjS[73] + tmpFx[93]*tmpObjS[86] + tmpFx[107]*tmpObjS[99] + tmpFx[121]*tmpObjS[112] + tmpFx[135]*tmpObjS[125] + tmpFx[149]*tmpObjS[138] + tmpFx[163]*tmpObjS[151] + tmpFx[177]*tmpObjS[164];
tmpQ2[126] = + tmpFx[9]*tmpObjS[9] + tmpFx[23]*tmpObjS[22] + tmpFx[37]*tmpObjS[35] + tmpFx[51]*tmpObjS[48] + tmpFx[65]*tmpObjS[61] + tmpFx[79]*tmpObjS[74] + tmpFx[93]*tmpObjS[87] + tmpFx[107]*tmpObjS[100] + tmpFx[121]*tmpObjS[113] + tmpFx[135]*tmpObjS[126] + tmpFx[149]*tmpObjS[139] + tmpFx[163]*tmpObjS[152] + tmpFx[177]*tmpObjS[165];
tmpQ2[127] = + tmpFx[9]*tmpObjS[10] + tmpFx[23]*tmpObjS[23] + tmpFx[37]*tmpObjS[36] + tmpFx[51]*tmpObjS[49] + tmpFx[65]*tmpObjS[62] + tmpFx[79]*tmpObjS[75] + tmpFx[93]*tmpObjS[88] + tmpFx[107]*tmpObjS[101] + tmpFx[121]*tmpObjS[114] + tmpFx[135]*tmpObjS[127] + tmpFx[149]*tmpObjS[140] + tmpFx[163]*tmpObjS[153] + tmpFx[177]*tmpObjS[166];
tmpQ2[128] = + tmpFx[9]*tmpObjS[11] + tmpFx[23]*tmpObjS[24] + tmpFx[37]*tmpObjS[37] + tmpFx[51]*tmpObjS[50] + tmpFx[65]*tmpObjS[63] + tmpFx[79]*tmpObjS[76] + tmpFx[93]*tmpObjS[89] + tmpFx[107]*tmpObjS[102] + tmpFx[121]*tmpObjS[115] + tmpFx[135]*tmpObjS[128] + tmpFx[149]*tmpObjS[141] + tmpFx[163]*tmpObjS[154] + tmpFx[177]*tmpObjS[167];
tmpQ2[129] = + tmpFx[9]*tmpObjS[12] + tmpFx[23]*tmpObjS[25] + tmpFx[37]*tmpObjS[38] + tmpFx[51]*tmpObjS[51] + tmpFx[65]*tmpObjS[64] + tmpFx[79]*tmpObjS[77] + tmpFx[93]*tmpObjS[90] + tmpFx[107]*tmpObjS[103] + tmpFx[121]*tmpObjS[116] + tmpFx[135]*tmpObjS[129] + tmpFx[149]*tmpObjS[142] + tmpFx[163]*tmpObjS[155] + tmpFx[177]*tmpObjS[168];
tmpQ2[130] = + tmpFx[10]*tmpObjS[0] + tmpFx[24]*tmpObjS[13] + tmpFx[38]*tmpObjS[26] + tmpFx[52]*tmpObjS[39] + tmpFx[66]*tmpObjS[52] + tmpFx[80]*tmpObjS[65] + tmpFx[94]*tmpObjS[78] + tmpFx[108]*tmpObjS[91] + tmpFx[122]*tmpObjS[104] + tmpFx[136]*tmpObjS[117] + tmpFx[150]*tmpObjS[130] + tmpFx[164]*tmpObjS[143] + tmpFx[178]*tmpObjS[156];
tmpQ2[131] = + tmpFx[10]*tmpObjS[1] + tmpFx[24]*tmpObjS[14] + tmpFx[38]*tmpObjS[27] + tmpFx[52]*tmpObjS[40] + tmpFx[66]*tmpObjS[53] + tmpFx[80]*tmpObjS[66] + tmpFx[94]*tmpObjS[79] + tmpFx[108]*tmpObjS[92] + tmpFx[122]*tmpObjS[105] + tmpFx[136]*tmpObjS[118] + tmpFx[150]*tmpObjS[131] + tmpFx[164]*tmpObjS[144] + tmpFx[178]*tmpObjS[157];
tmpQ2[132] = + tmpFx[10]*tmpObjS[2] + tmpFx[24]*tmpObjS[15] + tmpFx[38]*tmpObjS[28] + tmpFx[52]*tmpObjS[41] + tmpFx[66]*tmpObjS[54] + tmpFx[80]*tmpObjS[67] + tmpFx[94]*tmpObjS[80] + tmpFx[108]*tmpObjS[93] + tmpFx[122]*tmpObjS[106] + tmpFx[136]*tmpObjS[119] + tmpFx[150]*tmpObjS[132] + tmpFx[164]*tmpObjS[145] + tmpFx[178]*tmpObjS[158];
tmpQ2[133] = + tmpFx[10]*tmpObjS[3] + tmpFx[24]*tmpObjS[16] + tmpFx[38]*tmpObjS[29] + tmpFx[52]*tmpObjS[42] + tmpFx[66]*tmpObjS[55] + tmpFx[80]*tmpObjS[68] + tmpFx[94]*tmpObjS[81] + tmpFx[108]*tmpObjS[94] + tmpFx[122]*tmpObjS[107] + tmpFx[136]*tmpObjS[120] + tmpFx[150]*tmpObjS[133] + tmpFx[164]*tmpObjS[146] + tmpFx[178]*tmpObjS[159];
tmpQ2[134] = + tmpFx[10]*tmpObjS[4] + tmpFx[24]*tmpObjS[17] + tmpFx[38]*tmpObjS[30] + tmpFx[52]*tmpObjS[43] + tmpFx[66]*tmpObjS[56] + tmpFx[80]*tmpObjS[69] + tmpFx[94]*tmpObjS[82] + tmpFx[108]*tmpObjS[95] + tmpFx[122]*tmpObjS[108] + tmpFx[136]*tmpObjS[121] + tmpFx[150]*tmpObjS[134] + tmpFx[164]*tmpObjS[147] + tmpFx[178]*tmpObjS[160];
tmpQ2[135] = + tmpFx[10]*tmpObjS[5] + tmpFx[24]*tmpObjS[18] + tmpFx[38]*tmpObjS[31] + tmpFx[52]*tmpObjS[44] + tmpFx[66]*tmpObjS[57] + tmpFx[80]*tmpObjS[70] + tmpFx[94]*tmpObjS[83] + tmpFx[108]*tmpObjS[96] + tmpFx[122]*tmpObjS[109] + tmpFx[136]*tmpObjS[122] + tmpFx[150]*tmpObjS[135] + tmpFx[164]*tmpObjS[148] + tmpFx[178]*tmpObjS[161];
tmpQ2[136] = + tmpFx[10]*tmpObjS[6] + tmpFx[24]*tmpObjS[19] + tmpFx[38]*tmpObjS[32] + tmpFx[52]*tmpObjS[45] + tmpFx[66]*tmpObjS[58] + tmpFx[80]*tmpObjS[71] + tmpFx[94]*tmpObjS[84] + tmpFx[108]*tmpObjS[97] + tmpFx[122]*tmpObjS[110] + tmpFx[136]*tmpObjS[123] + tmpFx[150]*tmpObjS[136] + tmpFx[164]*tmpObjS[149] + tmpFx[178]*tmpObjS[162];
tmpQ2[137] = + tmpFx[10]*tmpObjS[7] + tmpFx[24]*tmpObjS[20] + tmpFx[38]*tmpObjS[33] + tmpFx[52]*tmpObjS[46] + tmpFx[66]*tmpObjS[59] + tmpFx[80]*tmpObjS[72] + tmpFx[94]*tmpObjS[85] + tmpFx[108]*tmpObjS[98] + tmpFx[122]*tmpObjS[111] + tmpFx[136]*tmpObjS[124] + tmpFx[150]*tmpObjS[137] + tmpFx[164]*tmpObjS[150] + tmpFx[178]*tmpObjS[163];
tmpQ2[138] = + tmpFx[10]*tmpObjS[8] + tmpFx[24]*tmpObjS[21] + tmpFx[38]*tmpObjS[34] + tmpFx[52]*tmpObjS[47] + tmpFx[66]*tmpObjS[60] + tmpFx[80]*tmpObjS[73] + tmpFx[94]*tmpObjS[86] + tmpFx[108]*tmpObjS[99] + tmpFx[122]*tmpObjS[112] + tmpFx[136]*tmpObjS[125] + tmpFx[150]*tmpObjS[138] + tmpFx[164]*tmpObjS[151] + tmpFx[178]*tmpObjS[164];
tmpQ2[139] = + tmpFx[10]*tmpObjS[9] + tmpFx[24]*tmpObjS[22] + tmpFx[38]*tmpObjS[35] + tmpFx[52]*tmpObjS[48] + tmpFx[66]*tmpObjS[61] + tmpFx[80]*tmpObjS[74] + tmpFx[94]*tmpObjS[87] + tmpFx[108]*tmpObjS[100] + tmpFx[122]*tmpObjS[113] + tmpFx[136]*tmpObjS[126] + tmpFx[150]*tmpObjS[139] + tmpFx[164]*tmpObjS[152] + tmpFx[178]*tmpObjS[165];
tmpQ2[140] = + tmpFx[10]*tmpObjS[10] + tmpFx[24]*tmpObjS[23] + tmpFx[38]*tmpObjS[36] + tmpFx[52]*tmpObjS[49] + tmpFx[66]*tmpObjS[62] + tmpFx[80]*tmpObjS[75] + tmpFx[94]*tmpObjS[88] + tmpFx[108]*tmpObjS[101] + tmpFx[122]*tmpObjS[114] + tmpFx[136]*tmpObjS[127] + tmpFx[150]*tmpObjS[140] + tmpFx[164]*tmpObjS[153] + tmpFx[178]*tmpObjS[166];
tmpQ2[141] = + tmpFx[10]*tmpObjS[11] + tmpFx[24]*tmpObjS[24] + tmpFx[38]*tmpObjS[37] + tmpFx[52]*tmpObjS[50] + tmpFx[66]*tmpObjS[63] + tmpFx[80]*tmpObjS[76] + tmpFx[94]*tmpObjS[89] + tmpFx[108]*tmpObjS[102] + tmpFx[122]*tmpObjS[115] + tmpFx[136]*tmpObjS[128] + tmpFx[150]*tmpObjS[141] + tmpFx[164]*tmpObjS[154] + tmpFx[178]*tmpObjS[167];
tmpQ2[142] = + tmpFx[10]*tmpObjS[12] + tmpFx[24]*tmpObjS[25] + tmpFx[38]*tmpObjS[38] + tmpFx[52]*tmpObjS[51] + tmpFx[66]*tmpObjS[64] + tmpFx[80]*tmpObjS[77] + tmpFx[94]*tmpObjS[90] + tmpFx[108]*tmpObjS[103] + tmpFx[122]*tmpObjS[116] + tmpFx[136]*tmpObjS[129] + tmpFx[150]*tmpObjS[142] + tmpFx[164]*tmpObjS[155] + tmpFx[178]*tmpObjS[168];
tmpQ2[143] = + tmpFx[11]*tmpObjS[0] + tmpFx[25]*tmpObjS[13] + tmpFx[39]*tmpObjS[26] + tmpFx[53]*tmpObjS[39] + tmpFx[67]*tmpObjS[52] + tmpFx[81]*tmpObjS[65] + tmpFx[95]*tmpObjS[78] + tmpFx[109]*tmpObjS[91] + tmpFx[123]*tmpObjS[104] + tmpFx[137]*tmpObjS[117] + tmpFx[151]*tmpObjS[130] + tmpFx[165]*tmpObjS[143] + tmpFx[179]*tmpObjS[156];
tmpQ2[144] = + tmpFx[11]*tmpObjS[1] + tmpFx[25]*tmpObjS[14] + tmpFx[39]*tmpObjS[27] + tmpFx[53]*tmpObjS[40] + tmpFx[67]*tmpObjS[53] + tmpFx[81]*tmpObjS[66] + tmpFx[95]*tmpObjS[79] + tmpFx[109]*tmpObjS[92] + tmpFx[123]*tmpObjS[105] + tmpFx[137]*tmpObjS[118] + tmpFx[151]*tmpObjS[131] + tmpFx[165]*tmpObjS[144] + tmpFx[179]*tmpObjS[157];
tmpQ2[145] = + tmpFx[11]*tmpObjS[2] + tmpFx[25]*tmpObjS[15] + tmpFx[39]*tmpObjS[28] + tmpFx[53]*tmpObjS[41] + tmpFx[67]*tmpObjS[54] + tmpFx[81]*tmpObjS[67] + tmpFx[95]*tmpObjS[80] + tmpFx[109]*tmpObjS[93] + tmpFx[123]*tmpObjS[106] + tmpFx[137]*tmpObjS[119] + tmpFx[151]*tmpObjS[132] + tmpFx[165]*tmpObjS[145] + tmpFx[179]*tmpObjS[158];
tmpQ2[146] = + tmpFx[11]*tmpObjS[3] + tmpFx[25]*tmpObjS[16] + tmpFx[39]*tmpObjS[29] + tmpFx[53]*tmpObjS[42] + tmpFx[67]*tmpObjS[55] + tmpFx[81]*tmpObjS[68] + tmpFx[95]*tmpObjS[81] + tmpFx[109]*tmpObjS[94] + tmpFx[123]*tmpObjS[107] + tmpFx[137]*tmpObjS[120] + tmpFx[151]*tmpObjS[133] + tmpFx[165]*tmpObjS[146] + tmpFx[179]*tmpObjS[159];
tmpQ2[147] = + tmpFx[11]*tmpObjS[4] + tmpFx[25]*tmpObjS[17] + tmpFx[39]*tmpObjS[30] + tmpFx[53]*tmpObjS[43] + tmpFx[67]*tmpObjS[56] + tmpFx[81]*tmpObjS[69] + tmpFx[95]*tmpObjS[82] + tmpFx[109]*tmpObjS[95] + tmpFx[123]*tmpObjS[108] + tmpFx[137]*tmpObjS[121] + tmpFx[151]*tmpObjS[134] + tmpFx[165]*tmpObjS[147] + tmpFx[179]*tmpObjS[160];
tmpQ2[148] = + tmpFx[11]*tmpObjS[5] + tmpFx[25]*tmpObjS[18] + tmpFx[39]*tmpObjS[31] + tmpFx[53]*tmpObjS[44] + tmpFx[67]*tmpObjS[57] + tmpFx[81]*tmpObjS[70] + tmpFx[95]*tmpObjS[83] + tmpFx[109]*tmpObjS[96] + tmpFx[123]*tmpObjS[109] + tmpFx[137]*tmpObjS[122] + tmpFx[151]*tmpObjS[135] + tmpFx[165]*tmpObjS[148] + tmpFx[179]*tmpObjS[161];
tmpQ2[149] = + tmpFx[11]*tmpObjS[6] + tmpFx[25]*tmpObjS[19] + tmpFx[39]*tmpObjS[32] + tmpFx[53]*tmpObjS[45] + tmpFx[67]*tmpObjS[58] + tmpFx[81]*tmpObjS[71] + tmpFx[95]*tmpObjS[84] + tmpFx[109]*tmpObjS[97] + tmpFx[123]*tmpObjS[110] + tmpFx[137]*tmpObjS[123] + tmpFx[151]*tmpObjS[136] + tmpFx[165]*tmpObjS[149] + tmpFx[179]*tmpObjS[162];
tmpQ2[150] = + tmpFx[11]*tmpObjS[7] + tmpFx[25]*tmpObjS[20] + tmpFx[39]*tmpObjS[33] + tmpFx[53]*tmpObjS[46] + tmpFx[67]*tmpObjS[59] + tmpFx[81]*tmpObjS[72] + tmpFx[95]*tmpObjS[85] + tmpFx[109]*tmpObjS[98] + tmpFx[123]*tmpObjS[111] + tmpFx[137]*tmpObjS[124] + tmpFx[151]*tmpObjS[137] + tmpFx[165]*tmpObjS[150] + tmpFx[179]*tmpObjS[163];
tmpQ2[151] = + tmpFx[11]*tmpObjS[8] + tmpFx[25]*tmpObjS[21] + tmpFx[39]*tmpObjS[34] + tmpFx[53]*tmpObjS[47] + tmpFx[67]*tmpObjS[60] + tmpFx[81]*tmpObjS[73] + tmpFx[95]*tmpObjS[86] + tmpFx[109]*tmpObjS[99] + tmpFx[123]*tmpObjS[112] + tmpFx[137]*tmpObjS[125] + tmpFx[151]*tmpObjS[138] + tmpFx[165]*tmpObjS[151] + tmpFx[179]*tmpObjS[164];
tmpQ2[152] = + tmpFx[11]*tmpObjS[9] + tmpFx[25]*tmpObjS[22] + tmpFx[39]*tmpObjS[35] + tmpFx[53]*tmpObjS[48] + tmpFx[67]*tmpObjS[61] + tmpFx[81]*tmpObjS[74] + tmpFx[95]*tmpObjS[87] + tmpFx[109]*tmpObjS[100] + tmpFx[123]*tmpObjS[113] + tmpFx[137]*tmpObjS[126] + tmpFx[151]*tmpObjS[139] + tmpFx[165]*tmpObjS[152] + tmpFx[179]*tmpObjS[165];
tmpQ2[153] = + tmpFx[11]*tmpObjS[10] + tmpFx[25]*tmpObjS[23] + tmpFx[39]*tmpObjS[36] + tmpFx[53]*tmpObjS[49] + tmpFx[67]*tmpObjS[62] + tmpFx[81]*tmpObjS[75] + tmpFx[95]*tmpObjS[88] + tmpFx[109]*tmpObjS[101] + tmpFx[123]*tmpObjS[114] + tmpFx[137]*tmpObjS[127] + tmpFx[151]*tmpObjS[140] + tmpFx[165]*tmpObjS[153] + tmpFx[179]*tmpObjS[166];
tmpQ2[154] = + tmpFx[11]*tmpObjS[11] + tmpFx[25]*tmpObjS[24] + tmpFx[39]*tmpObjS[37] + tmpFx[53]*tmpObjS[50] + tmpFx[67]*tmpObjS[63] + tmpFx[81]*tmpObjS[76] + tmpFx[95]*tmpObjS[89] + tmpFx[109]*tmpObjS[102] + tmpFx[123]*tmpObjS[115] + tmpFx[137]*tmpObjS[128] + tmpFx[151]*tmpObjS[141] + tmpFx[165]*tmpObjS[154] + tmpFx[179]*tmpObjS[167];
tmpQ2[155] = + tmpFx[11]*tmpObjS[12] + tmpFx[25]*tmpObjS[25] + tmpFx[39]*tmpObjS[38] + tmpFx[53]*tmpObjS[51] + tmpFx[67]*tmpObjS[64] + tmpFx[81]*tmpObjS[77] + tmpFx[95]*tmpObjS[90] + tmpFx[109]*tmpObjS[103] + tmpFx[123]*tmpObjS[116] + tmpFx[137]*tmpObjS[129] + tmpFx[151]*tmpObjS[142] + tmpFx[165]*tmpObjS[155] + tmpFx[179]*tmpObjS[168];
tmpQ2[156] = + tmpFx[12]*tmpObjS[0] + tmpFx[26]*tmpObjS[13] + tmpFx[40]*tmpObjS[26] + tmpFx[54]*tmpObjS[39] + tmpFx[68]*tmpObjS[52] + tmpFx[82]*tmpObjS[65] + tmpFx[96]*tmpObjS[78] + tmpFx[110]*tmpObjS[91] + tmpFx[124]*tmpObjS[104] + tmpFx[138]*tmpObjS[117] + tmpFx[152]*tmpObjS[130] + tmpFx[166]*tmpObjS[143] + tmpFx[180]*tmpObjS[156];
tmpQ2[157] = + tmpFx[12]*tmpObjS[1] + tmpFx[26]*tmpObjS[14] + tmpFx[40]*tmpObjS[27] + tmpFx[54]*tmpObjS[40] + tmpFx[68]*tmpObjS[53] + tmpFx[82]*tmpObjS[66] + tmpFx[96]*tmpObjS[79] + tmpFx[110]*tmpObjS[92] + tmpFx[124]*tmpObjS[105] + tmpFx[138]*tmpObjS[118] + tmpFx[152]*tmpObjS[131] + tmpFx[166]*tmpObjS[144] + tmpFx[180]*tmpObjS[157];
tmpQ2[158] = + tmpFx[12]*tmpObjS[2] + tmpFx[26]*tmpObjS[15] + tmpFx[40]*tmpObjS[28] + tmpFx[54]*tmpObjS[41] + tmpFx[68]*tmpObjS[54] + tmpFx[82]*tmpObjS[67] + tmpFx[96]*tmpObjS[80] + tmpFx[110]*tmpObjS[93] + tmpFx[124]*tmpObjS[106] + tmpFx[138]*tmpObjS[119] + tmpFx[152]*tmpObjS[132] + tmpFx[166]*tmpObjS[145] + tmpFx[180]*tmpObjS[158];
tmpQ2[159] = + tmpFx[12]*tmpObjS[3] + tmpFx[26]*tmpObjS[16] + tmpFx[40]*tmpObjS[29] + tmpFx[54]*tmpObjS[42] + tmpFx[68]*tmpObjS[55] + tmpFx[82]*tmpObjS[68] + tmpFx[96]*tmpObjS[81] + tmpFx[110]*tmpObjS[94] + tmpFx[124]*tmpObjS[107] + tmpFx[138]*tmpObjS[120] + tmpFx[152]*tmpObjS[133] + tmpFx[166]*tmpObjS[146] + tmpFx[180]*tmpObjS[159];
tmpQ2[160] = + tmpFx[12]*tmpObjS[4] + tmpFx[26]*tmpObjS[17] + tmpFx[40]*tmpObjS[30] + tmpFx[54]*tmpObjS[43] + tmpFx[68]*tmpObjS[56] + tmpFx[82]*tmpObjS[69] + tmpFx[96]*tmpObjS[82] + tmpFx[110]*tmpObjS[95] + tmpFx[124]*tmpObjS[108] + tmpFx[138]*tmpObjS[121] + tmpFx[152]*tmpObjS[134] + tmpFx[166]*tmpObjS[147] + tmpFx[180]*tmpObjS[160];
tmpQ2[161] = + tmpFx[12]*tmpObjS[5] + tmpFx[26]*tmpObjS[18] + tmpFx[40]*tmpObjS[31] + tmpFx[54]*tmpObjS[44] + tmpFx[68]*tmpObjS[57] + tmpFx[82]*tmpObjS[70] + tmpFx[96]*tmpObjS[83] + tmpFx[110]*tmpObjS[96] + tmpFx[124]*tmpObjS[109] + tmpFx[138]*tmpObjS[122] + tmpFx[152]*tmpObjS[135] + tmpFx[166]*tmpObjS[148] + tmpFx[180]*tmpObjS[161];
tmpQ2[162] = + tmpFx[12]*tmpObjS[6] + tmpFx[26]*tmpObjS[19] + tmpFx[40]*tmpObjS[32] + tmpFx[54]*tmpObjS[45] + tmpFx[68]*tmpObjS[58] + tmpFx[82]*tmpObjS[71] + tmpFx[96]*tmpObjS[84] + tmpFx[110]*tmpObjS[97] + tmpFx[124]*tmpObjS[110] + tmpFx[138]*tmpObjS[123] + tmpFx[152]*tmpObjS[136] + tmpFx[166]*tmpObjS[149] + tmpFx[180]*tmpObjS[162];
tmpQ2[163] = + tmpFx[12]*tmpObjS[7] + tmpFx[26]*tmpObjS[20] + tmpFx[40]*tmpObjS[33] + tmpFx[54]*tmpObjS[46] + tmpFx[68]*tmpObjS[59] + tmpFx[82]*tmpObjS[72] + tmpFx[96]*tmpObjS[85] + tmpFx[110]*tmpObjS[98] + tmpFx[124]*tmpObjS[111] + tmpFx[138]*tmpObjS[124] + tmpFx[152]*tmpObjS[137] + tmpFx[166]*tmpObjS[150] + tmpFx[180]*tmpObjS[163];
tmpQ2[164] = + tmpFx[12]*tmpObjS[8] + tmpFx[26]*tmpObjS[21] + tmpFx[40]*tmpObjS[34] + tmpFx[54]*tmpObjS[47] + tmpFx[68]*tmpObjS[60] + tmpFx[82]*tmpObjS[73] + tmpFx[96]*tmpObjS[86] + tmpFx[110]*tmpObjS[99] + tmpFx[124]*tmpObjS[112] + tmpFx[138]*tmpObjS[125] + tmpFx[152]*tmpObjS[138] + tmpFx[166]*tmpObjS[151] + tmpFx[180]*tmpObjS[164];
tmpQ2[165] = + tmpFx[12]*tmpObjS[9] + tmpFx[26]*tmpObjS[22] + tmpFx[40]*tmpObjS[35] + tmpFx[54]*tmpObjS[48] + tmpFx[68]*tmpObjS[61] + tmpFx[82]*tmpObjS[74] + tmpFx[96]*tmpObjS[87] + tmpFx[110]*tmpObjS[100] + tmpFx[124]*tmpObjS[113] + tmpFx[138]*tmpObjS[126] + tmpFx[152]*tmpObjS[139] + tmpFx[166]*tmpObjS[152] + tmpFx[180]*tmpObjS[165];
tmpQ2[166] = + tmpFx[12]*tmpObjS[10] + tmpFx[26]*tmpObjS[23] + tmpFx[40]*tmpObjS[36] + tmpFx[54]*tmpObjS[49] + tmpFx[68]*tmpObjS[62] + tmpFx[82]*tmpObjS[75] + tmpFx[96]*tmpObjS[88] + tmpFx[110]*tmpObjS[101] + tmpFx[124]*tmpObjS[114] + tmpFx[138]*tmpObjS[127] + tmpFx[152]*tmpObjS[140] + tmpFx[166]*tmpObjS[153] + tmpFx[180]*tmpObjS[166];
tmpQ2[167] = + tmpFx[12]*tmpObjS[11] + tmpFx[26]*tmpObjS[24] + tmpFx[40]*tmpObjS[37] + tmpFx[54]*tmpObjS[50] + tmpFx[68]*tmpObjS[63] + tmpFx[82]*tmpObjS[76] + tmpFx[96]*tmpObjS[89] + tmpFx[110]*tmpObjS[102] + tmpFx[124]*tmpObjS[115] + tmpFx[138]*tmpObjS[128] + tmpFx[152]*tmpObjS[141] + tmpFx[166]*tmpObjS[154] + tmpFx[180]*tmpObjS[167];
tmpQ2[168] = + tmpFx[12]*tmpObjS[12] + tmpFx[26]*tmpObjS[25] + tmpFx[40]*tmpObjS[38] + tmpFx[54]*tmpObjS[51] + tmpFx[68]*tmpObjS[64] + tmpFx[82]*tmpObjS[77] + tmpFx[96]*tmpObjS[90] + tmpFx[110]*tmpObjS[103] + tmpFx[124]*tmpObjS[116] + tmpFx[138]*tmpObjS[129] + tmpFx[152]*tmpObjS[142] + tmpFx[166]*tmpObjS[155] + tmpFx[180]*tmpObjS[168];
tmpQ2[169] = + tmpFx[13]*tmpObjS[0] + tmpFx[27]*tmpObjS[13] + tmpFx[41]*tmpObjS[26] + tmpFx[55]*tmpObjS[39] + tmpFx[69]*tmpObjS[52] + tmpFx[83]*tmpObjS[65] + tmpFx[97]*tmpObjS[78] + tmpFx[111]*tmpObjS[91] + tmpFx[125]*tmpObjS[104] + tmpFx[139]*tmpObjS[117] + tmpFx[153]*tmpObjS[130] + tmpFx[167]*tmpObjS[143] + tmpFx[181]*tmpObjS[156];
tmpQ2[170] = + tmpFx[13]*tmpObjS[1] + tmpFx[27]*tmpObjS[14] + tmpFx[41]*tmpObjS[27] + tmpFx[55]*tmpObjS[40] + tmpFx[69]*tmpObjS[53] + tmpFx[83]*tmpObjS[66] + tmpFx[97]*tmpObjS[79] + tmpFx[111]*tmpObjS[92] + tmpFx[125]*tmpObjS[105] + tmpFx[139]*tmpObjS[118] + tmpFx[153]*tmpObjS[131] + tmpFx[167]*tmpObjS[144] + tmpFx[181]*tmpObjS[157];
tmpQ2[171] = + tmpFx[13]*tmpObjS[2] + tmpFx[27]*tmpObjS[15] + tmpFx[41]*tmpObjS[28] + tmpFx[55]*tmpObjS[41] + tmpFx[69]*tmpObjS[54] + tmpFx[83]*tmpObjS[67] + tmpFx[97]*tmpObjS[80] + tmpFx[111]*tmpObjS[93] + tmpFx[125]*tmpObjS[106] + tmpFx[139]*tmpObjS[119] + tmpFx[153]*tmpObjS[132] + tmpFx[167]*tmpObjS[145] + tmpFx[181]*tmpObjS[158];
tmpQ2[172] = + tmpFx[13]*tmpObjS[3] + tmpFx[27]*tmpObjS[16] + tmpFx[41]*tmpObjS[29] + tmpFx[55]*tmpObjS[42] + tmpFx[69]*tmpObjS[55] + tmpFx[83]*tmpObjS[68] + tmpFx[97]*tmpObjS[81] + tmpFx[111]*tmpObjS[94] + tmpFx[125]*tmpObjS[107] + tmpFx[139]*tmpObjS[120] + tmpFx[153]*tmpObjS[133] + tmpFx[167]*tmpObjS[146] + tmpFx[181]*tmpObjS[159];
tmpQ2[173] = + tmpFx[13]*tmpObjS[4] + tmpFx[27]*tmpObjS[17] + tmpFx[41]*tmpObjS[30] + tmpFx[55]*tmpObjS[43] + tmpFx[69]*tmpObjS[56] + tmpFx[83]*tmpObjS[69] + tmpFx[97]*tmpObjS[82] + tmpFx[111]*tmpObjS[95] + tmpFx[125]*tmpObjS[108] + tmpFx[139]*tmpObjS[121] + tmpFx[153]*tmpObjS[134] + tmpFx[167]*tmpObjS[147] + tmpFx[181]*tmpObjS[160];
tmpQ2[174] = + tmpFx[13]*tmpObjS[5] + tmpFx[27]*tmpObjS[18] + tmpFx[41]*tmpObjS[31] + tmpFx[55]*tmpObjS[44] + tmpFx[69]*tmpObjS[57] + tmpFx[83]*tmpObjS[70] + tmpFx[97]*tmpObjS[83] + tmpFx[111]*tmpObjS[96] + tmpFx[125]*tmpObjS[109] + tmpFx[139]*tmpObjS[122] + tmpFx[153]*tmpObjS[135] + tmpFx[167]*tmpObjS[148] + tmpFx[181]*tmpObjS[161];
tmpQ2[175] = + tmpFx[13]*tmpObjS[6] + tmpFx[27]*tmpObjS[19] + tmpFx[41]*tmpObjS[32] + tmpFx[55]*tmpObjS[45] + tmpFx[69]*tmpObjS[58] + tmpFx[83]*tmpObjS[71] + tmpFx[97]*tmpObjS[84] + tmpFx[111]*tmpObjS[97] + tmpFx[125]*tmpObjS[110] + tmpFx[139]*tmpObjS[123] + tmpFx[153]*tmpObjS[136] + tmpFx[167]*tmpObjS[149] + tmpFx[181]*tmpObjS[162];
tmpQ2[176] = + tmpFx[13]*tmpObjS[7] + tmpFx[27]*tmpObjS[20] + tmpFx[41]*tmpObjS[33] + tmpFx[55]*tmpObjS[46] + tmpFx[69]*tmpObjS[59] + tmpFx[83]*tmpObjS[72] + tmpFx[97]*tmpObjS[85] + tmpFx[111]*tmpObjS[98] + tmpFx[125]*tmpObjS[111] + tmpFx[139]*tmpObjS[124] + tmpFx[153]*tmpObjS[137] + tmpFx[167]*tmpObjS[150] + tmpFx[181]*tmpObjS[163];
tmpQ2[177] = + tmpFx[13]*tmpObjS[8] + tmpFx[27]*tmpObjS[21] + tmpFx[41]*tmpObjS[34] + tmpFx[55]*tmpObjS[47] + tmpFx[69]*tmpObjS[60] + tmpFx[83]*tmpObjS[73] + tmpFx[97]*tmpObjS[86] + tmpFx[111]*tmpObjS[99] + tmpFx[125]*tmpObjS[112] + tmpFx[139]*tmpObjS[125] + tmpFx[153]*tmpObjS[138] + tmpFx[167]*tmpObjS[151] + tmpFx[181]*tmpObjS[164];
tmpQ2[178] = + tmpFx[13]*tmpObjS[9] + tmpFx[27]*tmpObjS[22] + tmpFx[41]*tmpObjS[35] + tmpFx[55]*tmpObjS[48] + tmpFx[69]*tmpObjS[61] + tmpFx[83]*tmpObjS[74] + tmpFx[97]*tmpObjS[87] + tmpFx[111]*tmpObjS[100] + tmpFx[125]*tmpObjS[113] + tmpFx[139]*tmpObjS[126] + tmpFx[153]*tmpObjS[139] + tmpFx[167]*tmpObjS[152] + tmpFx[181]*tmpObjS[165];
tmpQ2[179] = + tmpFx[13]*tmpObjS[10] + tmpFx[27]*tmpObjS[23] + tmpFx[41]*tmpObjS[36] + tmpFx[55]*tmpObjS[49] + tmpFx[69]*tmpObjS[62] + tmpFx[83]*tmpObjS[75] + tmpFx[97]*tmpObjS[88] + tmpFx[111]*tmpObjS[101] + tmpFx[125]*tmpObjS[114] + tmpFx[139]*tmpObjS[127] + tmpFx[153]*tmpObjS[140] + tmpFx[167]*tmpObjS[153] + tmpFx[181]*tmpObjS[166];
tmpQ2[180] = + tmpFx[13]*tmpObjS[11] + tmpFx[27]*tmpObjS[24] + tmpFx[41]*tmpObjS[37] + tmpFx[55]*tmpObjS[50] + tmpFx[69]*tmpObjS[63] + tmpFx[83]*tmpObjS[76] + tmpFx[97]*tmpObjS[89] + tmpFx[111]*tmpObjS[102] + tmpFx[125]*tmpObjS[115] + tmpFx[139]*tmpObjS[128] + tmpFx[153]*tmpObjS[141] + tmpFx[167]*tmpObjS[154] + tmpFx[181]*tmpObjS[167];
tmpQ2[181] = + tmpFx[13]*tmpObjS[12] + tmpFx[27]*tmpObjS[25] + tmpFx[41]*tmpObjS[38] + tmpFx[55]*tmpObjS[51] + tmpFx[69]*tmpObjS[64] + tmpFx[83]*tmpObjS[77] + tmpFx[97]*tmpObjS[90] + tmpFx[111]*tmpObjS[103] + tmpFx[125]*tmpObjS[116] + tmpFx[139]*tmpObjS[129] + tmpFx[153]*tmpObjS[142] + tmpFx[167]*tmpObjS[155] + tmpFx[181]*tmpObjS[168];
tmpQ1[0] = + tmpQ2[0]*tmpFx[0] + tmpQ2[1]*tmpFx[14] + tmpQ2[2]*tmpFx[28] + tmpQ2[3]*tmpFx[42] + tmpQ2[4]*tmpFx[56] + tmpQ2[5]*tmpFx[70] + tmpQ2[6]*tmpFx[84] + tmpQ2[7]*tmpFx[98] + tmpQ2[8]*tmpFx[112] + tmpQ2[9]*tmpFx[126] + tmpQ2[10]*tmpFx[140] + tmpQ2[11]*tmpFx[154] + tmpQ2[12]*tmpFx[168];
tmpQ1[1] = + tmpQ2[0]*tmpFx[1] + tmpQ2[1]*tmpFx[15] + tmpQ2[2]*tmpFx[29] + tmpQ2[3]*tmpFx[43] + tmpQ2[4]*tmpFx[57] + tmpQ2[5]*tmpFx[71] + tmpQ2[6]*tmpFx[85] + tmpQ2[7]*tmpFx[99] + tmpQ2[8]*tmpFx[113] + tmpQ2[9]*tmpFx[127] + tmpQ2[10]*tmpFx[141] + tmpQ2[11]*tmpFx[155] + tmpQ2[12]*tmpFx[169];
tmpQ1[2] = + tmpQ2[0]*tmpFx[2] + tmpQ2[1]*tmpFx[16] + tmpQ2[2]*tmpFx[30] + tmpQ2[3]*tmpFx[44] + tmpQ2[4]*tmpFx[58] + tmpQ2[5]*tmpFx[72] + tmpQ2[6]*tmpFx[86] + tmpQ2[7]*tmpFx[100] + tmpQ2[8]*tmpFx[114] + tmpQ2[9]*tmpFx[128] + tmpQ2[10]*tmpFx[142] + tmpQ2[11]*tmpFx[156] + tmpQ2[12]*tmpFx[170];
tmpQ1[3] = + tmpQ2[0]*tmpFx[3] + tmpQ2[1]*tmpFx[17] + tmpQ2[2]*tmpFx[31] + tmpQ2[3]*tmpFx[45] + tmpQ2[4]*tmpFx[59] + tmpQ2[5]*tmpFx[73] + tmpQ2[6]*tmpFx[87] + tmpQ2[7]*tmpFx[101] + tmpQ2[8]*tmpFx[115] + tmpQ2[9]*tmpFx[129] + tmpQ2[10]*tmpFx[143] + tmpQ2[11]*tmpFx[157] + tmpQ2[12]*tmpFx[171];
tmpQ1[4] = + tmpQ2[0]*tmpFx[4] + tmpQ2[1]*tmpFx[18] + tmpQ2[2]*tmpFx[32] + tmpQ2[3]*tmpFx[46] + tmpQ2[4]*tmpFx[60] + tmpQ2[5]*tmpFx[74] + tmpQ2[6]*tmpFx[88] + tmpQ2[7]*tmpFx[102] + tmpQ2[8]*tmpFx[116] + tmpQ2[9]*tmpFx[130] + tmpQ2[10]*tmpFx[144] + tmpQ2[11]*tmpFx[158] + tmpQ2[12]*tmpFx[172];
tmpQ1[5] = + tmpQ2[0]*tmpFx[5] + tmpQ2[1]*tmpFx[19] + tmpQ2[2]*tmpFx[33] + tmpQ2[3]*tmpFx[47] + tmpQ2[4]*tmpFx[61] + tmpQ2[5]*tmpFx[75] + tmpQ2[6]*tmpFx[89] + tmpQ2[7]*tmpFx[103] + tmpQ2[8]*tmpFx[117] + tmpQ2[9]*tmpFx[131] + tmpQ2[10]*tmpFx[145] + tmpQ2[11]*tmpFx[159] + tmpQ2[12]*tmpFx[173];
tmpQ1[6] = + tmpQ2[0]*tmpFx[6] + tmpQ2[1]*tmpFx[20] + tmpQ2[2]*tmpFx[34] + tmpQ2[3]*tmpFx[48] + tmpQ2[4]*tmpFx[62] + tmpQ2[5]*tmpFx[76] + tmpQ2[6]*tmpFx[90] + tmpQ2[7]*tmpFx[104] + tmpQ2[8]*tmpFx[118] + tmpQ2[9]*tmpFx[132] + tmpQ2[10]*tmpFx[146] + tmpQ2[11]*tmpFx[160] + tmpQ2[12]*tmpFx[174];
tmpQ1[7] = + tmpQ2[0]*tmpFx[7] + tmpQ2[1]*tmpFx[21] + tmpQ2[2]*tmpFx[35] + tmpQ2[3]*tmpFx[49] + tmpQ2[4]*tmpFx[63] + tmpQ2[5]*tmpFx[77] + tmpQ2[6]*tmpFx[91] + tmpQ2[7]*tmpFx[105] + tmpQ2[8]*tmpFx[119] + tmpQ2[9]*tmpFx[133] + tmpQ2[10]*tmpFx[147] + tmpQ2[11]*tmpFx[161] + tmpQ2[12]*tmpFx[175];
tmpQ1[8] = + tmpQ2[0]*tmpFx[8] + tmpQ2[1]*tmpFx[22] + tmpQ2[2]*tmpFx[36] + tmpQ2[3]*tmpFx[50] + tmpQ2[4]*tmpFx[64] + tmpQ2[5]*tmpFx[78] + tmpQ2[6]*tmpFx[92] + tmpQ2[7]*tmpFx[106] + tmpQ2[8]*tmpFx[120] + tmpQ2[9]*tmpFx[134] + tmpQ2[10]*tmpFx[148] + tmpQ2[11]*tmpFx[162] + tmpQ2[12]*tmpFx[176];
tmpQ1[9] = + tmpQ2[0]*tmpFx[9] + tmpQ2[1]*tmpFx[23] + tmpQ2[2]*tmpFx[37] + tmpQ2[3]*tmpFx[51] + tmpQ2[4]*tmpFx[65] + tmpQ2[5]*tmpFx[79] + tmpQ2[6]*tmpFx[93] + tmpQ2[7]*tmpFx[107] + tmpQ2[8]*tmpFx[121] + tmpQ2[9]*tmpFx[135] + tmpQ2[10]*tmpFx[149] + tmpQ2[11]*tmpFx[163] + tmpQ2[12]*tmpFx[177];
tmpQ1[10] = + tmpQ2[0]*tmpFx[10] + tmpQ2[1]*tmpFx[24] + tmpQ2[2]*tmpFx[38] + tmpQ2[3]*tmpFx[52] + tmpQ2[4]*tmpFx[66] + tmpQ2[5]*tmpFx[80] + tmpQ2[6]*tmpFx[94] + tmpQ2[7]*tmpFx[108] + tmpQ2[8]*tmpFx[122] + tmpQ2[9]*tmpFx[136] + tmpQ2[10]*tmpFx[150] + tmpQ2[11]*tmpFx[164] + tmpQ2[12]*tmpFx[178];
tmpQ1[11] = + tmpQ2[0]*tmpFx[11] + tmpQ2[1]*tmpFx[25] + tmpQ2[2]*tmpFx[39] + tmpQ2[3]*tmpFx[53] + tmpQ2[4]*tmpFx[67] + tmpQ2[5]*tmpFx[81] + tmpQ2[6]*tmpFx[95] + tmpQ2[7]*tmpFx[109] + tmpQ2[8]*tmpFx[123] + tmpQ2[9]*tmpFx[137] + tmpQ2[10]*tmpFx[151] + tmpQ2[11]*tmpFx[165] + tmpQ2[12]*tmpFx[179];
tmpQ1[12] = + tmpQ2[0]*tmpFx[12] + tmpQ2[1]*tmpFx[26] + tmpQ2[2]*tmpFx[40] + tmpQ2[3]*tmpFx[54] + tmpQ2[4]*tmpFx[68] + tmpQ2[5]*tmpFx[82] + tmpQ2[6]*tmpFx[96] + tmpQ2[7]*tmpFx[110] + tmpQ2[8]*tmpFx[124] + tmpQ2[9]*tmpFx[138] + tmpQ2[10]*tmpFx[152] + tmpQ2[11]*tmpFx[166] + tmpQ2[12]*tmpFx[180];
tmpQ1[13] = + tmpQ2[0]*tmpFx[13] + tmpQ2[1]*tmpFx[27] + tmpQ2[2]*tmpFx[41] + tmpQ2[3]*tmpFx[55] + tmpQ2[4]*tmpFx[69] + tmpQ2[5]*tmpFx[83] + tmpQ2[6]*tmpFx[97] + tmpQ2[7]*tmpFx[111] + tmpQ2[8]*tmpFx[125] + tmpQ2[9]*tmpFx[139] + tmpQ2[10]*tmpFx[153] + tmpQ2[11]*tmpFx[167] + tmpQ2[12]*tmpFx[181];
tmpQ1[14] = + tmpQ2[13]*tmpFx[0] + tmpQ2[14]*tmpFx[14] + tmpQ2[15]*tmpFx[28] + tmpQ2[16]*tmpFx[42] + tmpQ2[17]*tmpFx[56] + tmpQ2[18]*tmpFx[70] + tmpQ2[19]*tmpFx[84] + tmpQ2[20]*tmpFx[98] + tmpQ2[21]*tmpFx[112] + tmpQ2[22]*tmpFx[126] + tmpQ2[23]*tmpFx[140] + tmpQ2[24]*tmpFx[154] + tmpQ2[25]*tmpFx[168];
tmpQ1[15] = + tmpQ2[13]*tmpFx[1] + tmpQ2[14]*tmpFx[15] + tmpQ2[15]*tmpFx[29] + tmpQ2[16]*tmpFx[43] + tmpQ2[17]*tmpFx[57] + tmpQ2[18]*tmpFx[71] + tmpQ2[19]*tmpFx[85] + tmpQ2[20]*tmpFx[99] + tmpQ2[21]*tmpFx[113] + tmpQ2[22]*tmpFx[127] + tmpQ2[23]*tmpFx[141] + tmpQ2[24]*tmpFx[155] + tmpQ2[25]*tmpFx[169];
tmpQ1[16] = + tmpQ2[13]*tmpFx[2] + tmpQ2[14]*tmpFx[16] + tmpQ2[15]*tmpFx[30] + tmpQ2[16]*tmpFx[44] + tmpQ2[17]*tmpFx[58] + tmpQ2[18]*tmpFx[72] + tmpQ2[19]*tmpFx[86] + tmpQ2[20]*tmpFx[100] + tmpQ2[21]*tmpFx[114] + tmpQ2[22]*tmpFx[128] + tmpQ2[23]*tmpFx[142] + tmpQ2[24]*tmpFx[156] + tmpQ2[25]*tmpFx[170];
tmpQ1[17] = + tmpQ2[13]*tmpFx[3] + tmpQ2[14]*tmpFx[17] + tmpQ2[15]*tmpFx[31] + tmpQ2[16]*tmpFx[45] + tmpQ2[17]*tmpFx[59] + tmpQ2[18]*tmpFx[73] + tmpQ2[19]*tmpFx[87] + tmpQ2[20]*tmpFx[101] + tmpQ2[21]*tmpFx[115] + tmpQ2[22]*tmpFx[129] + tmpQ2[23]*tmpFx[143] + tmpQ2[24]*tmpFx[157] + tmpQ2[25]*tmpFx[171];
tmpQ1[18] = + tmpQ2[13]*tmpFx[4] + tmpQ2[14]*tmpFx[18] + tmpQ2[15]*tmpFx[32] + tmpQ2[16]*tmpFx[46] + tmpQ2[17]*tmpFx[60] + tmpQ2[18]*tmpFx[74] + tmpQ2[19]*tmpFx[88] + tmpQ2[20]*tmpFx[102] + tmpQ2[21]*tmpFx[116] + tmpQ2[22]*tmpFx[130] + tmpQ2[23]*tmpFx[144] + tmpQ2[24]*tmpFx[158] + tmpQ2[25]*tmpFx[172];
tmpQ1[19] = + tmpQ2[13]*tmpFx[5] + tmpQ2[14]*tmpFx[19] + tmpQ2[15]*tmpFx[33] + tmpQ2[16]*tmpFx[47] + tmpQ2[17]*tmpFx[61] + tmpQ2[18]*tmpFx[75] + tmpQ2[19]*tmpFx[89] + tmpQ2[20]*tmpFx[103] + tmpQ2[21]*tmpFx[117] + tmpQ2[22]*tmpFx[131] + tmpQ2[23]*tmpFx[145] + tmpQ2[24]*tmpFx[159] + tmpQ2[25]*tmpFx[173];
tmpQ1[20] = + tmpQ2[13]*tmpFx[6] + tmpQ2[14]*tmpFx[20] + tmpQ2[15]*tmpFx[34] + tmpQ2[16]*tmpFx[48] + tmpQ2[17]*tmpFx[62] + tmpQ2[18]*tmpFx[76] + tmpQ2[19]*tmpFx[90] + tmpQ2[20]*tmpFx[104] + tmpQ2[21]*tmpFx[118] + tmpQ2[22]*tmpFx[132] + tmpQ2[23]*tmpFx[146] + tmpQ2[24]*tmpFx[160] + tmpQ2[25]*tmpFx[174];
tmpQ1[21] = + tmpQ2[13]*tmpFx[7] + tmpQ2[14]*tmpFx[21] + tmpQ2[15]*tmpFx[35] + tmpQ2[16]*tmpFx[49] + tmpQ2[17]*tmpFx[63] + tmpQ2[18]*tmpFx[77] + tmpQ2[19]*tmpFx[91] + tmpQ2[20]*tmpFx[105] + tmpQ2[21]*tmpFx[119] + tmpQ2[22]*tmpFx[133] + tmpQ2[23]*tmpFx[147] + tmpQ2[24]*tmpFx[161] + tmpQ2[25]*tmpFx[175];
tmpQ1[22] = + tmpQ2[13]*tmpFx[8] + tmpQ2[14]*tmpFx[22] + tmpQ2[15]*tmpFx[36] + tmpQ2[16]*tmpFx[50] + tmpQ2[17]*tmpFx[64] + tmpQ2[18]*tmpFx[78] + tmpQ2[19]*tmpFx[92] + tmpQ2[20]*tmpFx[106] + tmpQ2[21]*tmpFx[120] + tmpQ2[22]*tmpFx[134] + tmpQ2[23]*tmpFx[148] + tmpQ2[24]*tmpFx[162] + tmpQ2[25]*tmpFx[176];
tmpQ1[23] = + tmpQ2[13]*tmpFx[9] + tmpQ2[14]*tmpFx[23] + tmpQ2[15]*tmpFx[37] + tmpQ2[16]*tmpFx[51] + tmpQ2[17]*tmpFx[65] + tmpQ2[18]*tmpFx[79] + tmpQ2[19]*tmpFx[93] + tmpQ2[20]*tmpFx[107] + tmpQ2[21]*tmpFx[121] + tmpQ2[22]*tmpFx[135] + tmpQ2[23]*tmpFx[149] + tmpQ2[24]*tmpFx[163] + tmpQ2[25]*tmpFx[177];
tmpQ1[24] = + tmpQ2[13]*tmpFx[10] + tmpQ2[14]*tmpFx[24] + tmpQ2[15]*tmpFx[38] + tmpQ2[16]*tmpFx[52] + tmpQ2[17]*tmpFx[66] + tmpQ2[18]*tmpFx[80] + tmpQ2[19]*tmpFx[94] + tmpQ2[20]*tmpFx[108] + tmpQ2[21]*tmpFx[122] + tmpQ2[22]*tmpFx[136] + tmpQ2[23]*tmpFx[150] + tmpQ2[24]*tmpFx[164] + tmpQ2[25]*tmpFx[178];
tmpQ1[25] = + tmpQ2[13]*tmpFx[11] + tmpQ2[14]*tmpFx[25] + tmpQ2[15]*tmpFx[39] + tmpQ2[16]*tmpFx[53] + tmpQ2[17]*tmpFx[67] + tmpQ2[18]*tmpFx[81] + tmpQ2[19]*tmpFx[95] + tmpQ2[20]*tmpFx[109] + tmpQ2[21]*tmpFx[123] + tmpQ2[22]*tmpFx[137] + tmpQ2[23]*tmpFx[151] + tmpQ2[24]*tmpFx[165] + tmpQ2[25]*tmpFx[179];
tmpQ1[26] = + tmpQ2[13]*tmpFx[12] + tmpQ2[14]*tmpFx[26] + tmpQ2[15]*tmpFx[40] + tmpQ2[16]*tmpFx[54] + tmpQ2[17]*tmpFx[68] + tmpQ2[18]*tmpFx[82] + tmpQ2[19]*tmpFx[96] + tmpQ2[20]*tmpFx[110] + tmpQ2[21]*tmpFx[124] + tmpQ2[22]*tmpFx[138] + tmpQ2[23]*tmpFx[152] + tmpQ2[24]*tmpFx[166] + tmpQ2[25]*tmpFx[180];
tmpQ1[27] = + tmpQ2[13]*tmpFx[13] + tmpQ2[14]*tmpFx[27] + tmpQ2[15]*tmpFx[41] + tmpQ2[16]*tmpFx[55] + tmpQ2[17]*tmpFx[69] + tmpQ2[18]*tmpFx[83] + tmpQ2[19]*tmpFx[97] + tmpQ2[20]*tmpFx[111] + tmpQ2[21]*tmpFx[125] + tmpQ2[22]*tmpFx[139] + tmpQ2[23]*tmpFx[153] + tmpQ2[24]*tmpFx[167] + tmpQ2[25]*tmpFx[181];
tmpQ1[28] = + tmpQ2[26]*tmpFx[0] + tmpQ2[27]*tmpFx[14] + tmpQ2[28]*tmpFx[28] + tmpQ2[29]*tmpFx[42] + tmpQ2[30]*tmpFx[56] + tmpQ2[31]*tmpFx[70] + tmpQ2[32]*tmpFx[84] + tmpQ2[33]*tmpFx[98] + tmpQ2[34]*tmpFx[112] + tmpQ2[35]*tmpFx[126] + tmpQ2[36]*tmpFx[140] + tmpQ2[37]*tmpFx[154] + tmpQ2[38]*tmpFx[168];
tmpQ1[29] = + tmpQ2[26]*tmpFx[1] + tmpQ2[27]*tmpFx[15] + tmpQ2[28]*tmpFx[29] + tmpQ2[29]*tmpFx[43] + tmpQ2[30]*tmpFx[57] + tmpQ2[31]*tmpFx[71] + tmpQ2[32]*tmpFx[85] + tmpQ2[33]*tmpFx[99] + tmpQ2[34]*tmpFx[113] + tmpQ2[35]*tmpFx[127] + tmpQ2[36]*tmpFx[141] + tmpQ2[37]*tmpFx[155] + tmpQ2[38]*tmpFx[169];
tmpQ1[30] = + tmpQ2[26]*tmpFx[2] + tmpQ2[27]*tmpFx[16] + tmpQ2[28]*tmpFx[30] + tmpQ2[29]*tmpFx[44] + tmpQ2[30]*tmpFx[58] + tmpQ2[31]*tmpFx[72] + tmpQ2[32]*tmpFx[86] + tmpQ2[33]*tmpFx[100] + tmpQ2[34]*tmpFx[114] + tmpQ2[35]*tmpFx[128] + tmpQ2[36]*tmpFx[142] + tmpQ2[37]*tmpFx[156] + tmpQ2[38]*tmpFx[170];
tmpQ1[31] = + tmpQ2[26]*tmpFx[3] + tmpQ2[27]*tmpFx[17] + tmpQ2[28]*tmpFx[31] + tmpQ2[29]*tmpFx[45] + tmpQ2[30]*tmpFx[59] + tmpQ2[31]*tmpFx[73] + tmpQ2[32]*tmpFx[87] + tmpQ2[33]*tmpFx[101] + tmpQ2[34]*tmpFx[115] + tmpQ2[35]*tmpFx[129] + tmpQ2[36]*tmpFx[143] + tmpQ2[37]*tmpFx[157] + tmpQ2[38]*tmpFx[171];
tmpQ1[32] = + tmpQ2[26]*tmpFx[4] + tmpQ2[27]*tmpFx[18] + tmpQ2[28]*tmpFx[32] + tmpQ2[29]*tmpFx[46] + tmpQ2[30]*tmpFx[60] + tmpQ2[31]*tmpFx[74] + tmpQ2[32]*tmpFx[88] + tmpQ2[33]*tmpFx[102] + tmpQ2[34]*tmpFx[116] + tmpQ2[35]*tmpFx[130] + tmpQ2[36]*tmpFx[144] + tmpQ2[37]*tmpFx[158] + tmpQ2[38]*tmpFx[172];
tmpQ1[33] = + tmpQ2[26]*tmpFx[5] + tmpQ2[27]*tmpFx[19] + tmpQ2[28]*tmpFx[33] + tmpQ2[29]*tmpFx[47] + tmpQ2[30]*tmpFx[61] + tmpQ2[31]*tmpFx[75] + tmpQ2[32]*tmpFx[89] + tmpQ2[33]*tmpFx[103] + tmpQ2[34]*tmpFx[117] + tmpQ2[35]*tmpFx[131] + tmpQ2[36]*tmpFx[145] + tmpQ2[37]*tmpFx[159] + tmpQ2[38]*tmpFx[173];
tmpQ1[34] = + tmpQ2[26]*tmpFx[6] + tmpQ2[27]*tmpFx[20] + tmpQ2[28]*tmpFx[34] + tmpQ2[29]*tmpFx[48] + tmpQ2[30]*tmpFx[62] + tmpQ2[31]*tmpFx[76] + tmpQ2[32]*tmpFx[90] + tmpQ2[33]*tmpFx[104] + tmpQ2[34]*tmpFx[118] + tmpQ2[35]*tmpFx[132] + tmpQ2[36]*tmpFx[146] + tmpQ2[37]*tmpFx[160] + tmpQ2[38]*tmpFx[174];
tmpQ1[35] = + tmpQ2[26]*tmpFx[7] + tmpQ2[27]*tmpFx[21] + tmpQ2[28]*tmpFx[35] + tmpQ2[29]*tmpFx[49] + tmpQ2[30]*tmpFx[63] + tmpQ2[31]*tmpFx[77] + tmpQ2[32]*tmpFx[91] + tmpQ2[33]*tmpFx[105] + tmpQ2[34]*tmpFx[119] + tmpQ2[35]*tmpFx[133] + tmpQ2[36]*tmpFx[147] + tmpQ2[37]*tmpFx[161] + tmpQ2[38]*tmpFx[175];
tmpQ1[36] = + tmpQ2[26]*tmpFx[8] + tmpQ2[27]*tmpFx[22] + tmpQ2[28]*tmpFx[36] + tmpQ2[29]*tmpFx[50] + tmpQ2[30]*tmpFx[64] + tmpQ2[31]*tmpFx[78] + tmpQ2[32]*tmpFx[92] + tmpQ2[33]*tmpFx[106] + tmpQ2[34]*tmpFx[120] + tmpQ2[35]*tmpFx[134] + tmpQ2[36]*tmpFx[148] + tmpQ2[37]*tmpFx[162] + tmpQ2[38]*tmpFx[176];
tmpQ1[37] = + tmpQ2[26]*tmpFx[9] + tmpQ2[27]*tmpFx[23] + tmpQ2[28]*tmpFx[37] + tmpQ2[29]*tmpFx[51] + tmpQ2[30]*tmpFx[65] + tmpQ2[31]*tmpFx[79] + tmpQ2[32]*tmpFx[93] + tmpQ2[33]*tmpFx[107] + tmpQ2[34]*tmpFx[121] + tmpQ2[35]*tmpFx[135] + tmpQ2[36]*tmpFx[149] + tmpQ2[37]*tmpFx[163] + tmpQ2[38]*tmpFx[177];
tmpQ1[38] = + tmpQ2[26]*tmpFx[10] + tmpQ2[27]*tmpFx[24] + tmpQ2[28]*tmpFx[38] + tmpQ2[29]*tmpFx[52] + tmpQ2[30]*tmpFx[66] + tmpQ2[31]*tmpFx[80] + tmpQ2[32]*tmpFx[94] + tmpQ2[33]*tmpFx[108] + tmpQ2[34]*tmpFx[122] + tmpQ2[35]*tmpFx[136] + tmpQ2[36]*tmpFx[150] + tmpQ2[37]*tmpFx[164] + tmpQ2[38]*tmpFx[178];
tmpQ1[39] = + tmpQ2[26]*tmpFx[11] + tmpQ2[27]*tmpFx[25] + tmpQ2[28]*tmpFx[39] + tmpQ2[29]*tmpFx[53] + tmpQ2[30]*tmpFx[67] + tmpQ2[31]*tmpFx[81] + tmpQ2[32]*tmpFx[95] + tmpQ2[33]*tmpFx[109] + tmpQ2[34]*tmpFx[123] + tmpQ2[35]*tmpFx[137] + tmpQ2[36]*tmpFx[151] + tmpQ2[37]*tmpFx[165] + tmpQ2[38]*tmpFx[179];
tmpQ1[40] = + tmpQ2[26]*tmpFx[12] + tmpQ2[27]*tmpFx[26] + tmpQ2[28]*tmpFx[40] + tmpQ2[29]*tmpFx[54] + tmpQ2[30]*tmpFx[68] + tmpQ2[31]*tmpFx[82] + tmpQ2[32]*tmpFx[96] + tmpQ2[33]*tmpFx[110] + tmpQ2[34]*tmpFx[124] + tmpQ2[35]*tmpFx[138] + tmpQ2[36]*tmpFx[152] + tmpQ2[37]*tmpFx[166] + tmpQ2[38]*tmpFx[180];
tmpQ1[41] = + tmpQ2[26]*tmpFx[13] + tmpQ2[27]*tmpFx[27] + tmpQ2[28]*tmpFx[41] + tmpQ2[29]*tmpFx[55] + tmpQ2[30]*tmpFx[69] + tmpQ2[31]*tmpFx[83] + tmpQ2[32]*tmpFx[97] + tmpQ2[33]*tmpFx[111] + tmpQ2[34]*tmpFx[125] + tmpQ2[35]*tmpFx[139] + tmpQ2[36]*tmpFx[153] + tmpQ2[37]*tmpFx[167] + tmpQ2[38]*tmpFx[181];
tmpQ1[42] = + tmpQ2[39]*tmpFx[0] + tmpQ2[40]*tmpFx[14] + tmpQ2[41]*tmpFx[28] + tmpQ2[42]*tmpFx[42] + tmpQ2[43]*tmpFx[56] + tmpQ2[44]*tmpFx[70] + tmpQ2[45]*tmpFx[84] + tmpQ2[46]*tmpFx[98] + tmpQ2[47]*tmpFx[112] + tmpQ2[48]*tmpFx[126] + tmpQ2[49]*tmpFx[140] + tmpQ2[50]*tmpFx[154] + tmpQ2[51]*tmpFx[168];
tmpQ1[43] = + tmpQ2[39]*tmpFx[1] + tmpQ2[40]*tmpFx[15] + tmpQ2[41]*tmpFx[29] + tmpQ2[42]*tmpFx[43] + tmpQ2[43]*tmpFx[57] + tmpQ2[44]*tmpFx[71] + tmpQ2[45]*tmpFx[85] + tmpQ2[46]*tmpFx[99] + tmpQ2[47]*tmpFx[113] + tmpQ2[48]*tmpFx[127] + tmpQ2[49]*tmpFx[141] + tmpQ2[50]*tmpFx[155] + tmpQ2[51]*tmpFx[169];
tmpQ1[44] = + tmpQ2[39]*tmpFx[2] + tmpQ2[40]*tmpFx[16] + tmpQ2[41]*tmpFx[30] + tmpQ2[42]*tmpFx[44] + tmpQ2[43]*tmpFx[58] + tmpQ2[44]*tmpFx[72] + tmpQ2[45]*tmpFx[86] + tmpQ2[46]*tmpFx[100] + tmpQ2[47]*tmpFx[114] + tmpQ2[48]*tmpFx[128] + tmpQ2[49]*tmpFx[142] + tmpQ2[50]*tmpFx[156] + tmpQ2[51]*tmpFx[170];
tmpQ1[45] = + tmpQ2[39]*tmpFx[3] + tmpQ2[40]*tmpFx[17] + tmpQ2[41]*tmpFx[31] + tmpQ2[42]*tmpFx[45] + tmpQ2[43]*tmpFx[59] + tmpQ2[44]*tmpFx[73] + tmpQ2[45]*tmpFx[87] + tmpQ2[46]*tmpFx[101] + tmpQ2[47]*tmpFx[115] + tmpQ2[48]*tmpFx[129] + tmpQ2[49]*tmpFx[143] + tmpQ2[50]*tmpFx[157] + tmpQ2[51]*tmpFx[171];
tmpQ1[46] = + tmpQ2[39]*tmpFx[4] + tmpQ2[40]*tmpFx[18] + tmpQ2[41]*tmpFx[32] + tmpQ2[42]*tmpFx[46] + tmpQ2[43]*tmpFx[60] + tmpQ2[44]*tmpFx[74] + tmpQ2[45]*tmpFx[88] + tmpQ2[46]*tmpFx[102] + tmpQ2[47]*tmpFx[116] + tmpQ2[48]*tmpFx[130] + tmpQ2[49]*tmpFx[144] + tmpQ2[50]*tmpFx[158] + tmpQ2[51]*tmpFx[172];
tmpQ1[47] = + tmpQ2[39]*tmpFx[5] + tmpQ2[40]*tmpFx[19] + tmpQ2[41]*tmpFx[33] + tmpQ2[42]*tmpFx[47] + tmpQ2[43]*tmpFx[61] + tmpQ2[44]*tmpFx[75] + tmpQ2[45]*tmpFx[89] + tmpQ2[46]*tmpFx[103] + tmpQ2[47]*tmpFx[117] + tmpQ2[48]*tmpFx[131] + tmpQ2[49]*tmpFx[145] + tmpQ2[50]*tmpFx[159] + tmpQ2[51]*tmpFx[173];
tmpQ1[48] = + tmpQ2[39]*tmpFx[6] + tmpQ2[40]*tmpFx[20] + tmpQ2[41]*tmpFx[34] + tmpQ2[42]*tmpFx[48] + tmpQ2[43]*tmpFx[62] + tmpQ2[44]*tmpFx[76] + tmpQ2[45]*tmpFx[90] + tmpQ2[46]*tmpFx[104] + tmpQ2[47]*tmpFx[118] + tmpQ2[48]*tmpFx[132] + tmpQ2[49]*tmpFx[146] + tmpQ2[50]*tmpFx[160] + tmpQ2[51]*tmpFx[174];
tmpQ1[49] = + tmpQ2[39]*tmpFx[7] + tmpQ2[40]*tmpFx[21] + tmpQ2[41]*tmpFx[35] + tmpQ2[42]*tmpFx[49] + tmpQ2[43]*tmpFx[63] + tmpQ2[44]*tmpFx[77] + tmpQ2[45]*tmpFx[91] + tmpQ2[46]*tmpFx[105] + tmpQ2[47]*tmpFx[119] + tmpQ2[48]*tmpFx[133] + tmpQ2[49]*tmpFx[147] + tmpQ2[50]*tmpFx[161] + tmpQ2[51]*tmpFx[175];
tmpQ1[50] = + tmpQ2[39]*tmpFx[8] + tmpQ2[40]*tmpFx[22] + tmpQ2[41]*tmpFx[36] + tmpQ2[42]*tmpFx[50] + tmpQ2[43]*tmpFx[64] + tmpQ2[44]*tmpFx[78] + tmpQ2[45]*tmpFx[92] + tmpQ2[46]*tmpFx[106] + tmpQ2[47]*tmpFx[120] + tmpQ2[48]*tmpFx[134] + tmpQ2[49]*tmpFx[148] + tmpQ2[50]*tmpFx[162] + tmpQ2[51]*tmpFx[176];
tmpQ1[51] = + tmpQ2[39]*tmpFx[9] + tmpQ2[40]*tmpFx[23] + tmpQ2[41]*tmpFx[37] + tmpQ2[42]*tmpFx[51] + tmpQ2[43]*tmpFx[65] + tmpQ2[44]*tmpFx[79] + tmpQ2[45]*tmpFx[93] + tmpQ2[46]*tmpFx[107] + tmpQ2[47]*tmpFx[121] + tmpQ2[48]*tmpFx[135] + tmpQ2[49]*tmpFx[149] + tmpQ2[50]*tmpFx[163] + tmpQ2[51]*tmpFx[177];
tmpQ1[52] = + tmpQ2[39]*tmpFx[10] + tmpQ2[40]*tmpFx[24] + tmpQ2[41]*tmpFx[38] + tmpQ2[42]*tmpFx[52] + tmpQ2[43]*tmpFx[66] + tmpQ2[44]*tmpFx[80] + tmpQ2[45]*tmpFx[94] + tmpQ2[46]*tmpFx[108] + tmpQ2[47]*tmpFx[122] + tmpQ2[48]*tmpFx[136] + tmpQ2[49]*tmpFx[150] + tmpQ2[50]*tmpFx[164] + tmpQ2[51]*tmpFx[178];
tmpQ1[53] = + tmpQ2[39]*tmpFx[11] + tmpQ2[40]*tmpFx[25] + tmpQ2[41]*tmpFx[39] + tmpQ2[42]*tmpFx[53] + tmpQ2[43]*tmpFx[67] + tmpQ2[44]*tmpFx[81] + tmpQ2[45]*tmpFx[95] + tmpQ2[46]*tmpFx[109] + tmpQ2[47]*tmpFx[123] + tmpQ2[48]*tmpFx[137] + tmpQ2[49]*tmpFx[151] + tmpQ2[50]*tmpFx[165] + tmpQ2[51]*tmpFx[179];
tmpQ1[54] = + tmpQ2[39]*tmpFx[12] + tmpQ2[40]*tmpFx[26] + tmpQ2[41]*tmpFx[40] + tmpQ2[42]*tmpFx[54] + tmpQ2[43]*tmpFx[68] + tmpQ2[44]*tmpFx[82] + tmpQ2[45]*tmpFx[96] + tmpQ2[46]*tmpFx[110] + tmpQ2[47]*tmpFx[124] + tmpQ2[48]*tmpFx[138] + tmpQ2[49]*tmpFx[152] + tmpQ2[50]*tmpFx[166] + tmpQ2[51]*tmpFx[180];
tmpQ1[55] = + tmpQ2[39]*tmpFx[13] + tmpQ2[40]*tmpFx[27] + tmpQ2[41]*tmpFx[41] + tmpQ2[42]*tmpFx[55] + tmpQ2[43]*tmpFx[69] + tmpQ2[44]*tmpFx[83] + tmpQ2[45]*tmpFx[97] + tmpQ2[46]*tmpFx[111] + tmpQ2[47]*tmpFx[125] + tmpQ2[48]*tmpFx[139] + tmpQ2[49]*tmpFx[153] + tmpQ2[50]*tmpFx[167] + tmpQ2[51]*tmpFx[181];
tmpQ1[56] = + tmpQ2[52]*tmpFx[0] + tmpQ2[53]*tmpFx[14] + tmpQ2[54]*tmpFx[28] + tmpQ2[55]*tmpFx[42] + tmpQ2[56]*tmpFx[56] + tmpQ2[57]*tmpFx[70] + tmpQ2[58]*tmpFx[84] + tmpQ2[59]*tmpFx[98] + tmpQ2[60]*tmpFx[112] + tmpQ2[61]*tmpFx[126] + tmpQ2[62]*tmpFx[140] + tmpQ2[63]*tmpFx[154] + tmpQ2[64]*tmpFx[168];
tmpQ1[57] = + tmpQ2[52]*tmpFx[1] + tmpQ2[53]*tmpFx[15] + tmpQ2[54]*tmpFx[29] + tmpQ2[55]*tmpFx[43] + tmpQ2[56]*tmpFx[57] + tmpQ2[57]*tmpFx[71] + tmpQ2[58]*tmpFx[85] + tmpQ2[59]*tmpFx[99] + tmpQ2[60]*tmpFx[113] + tmpQ2[61]*tmpFx[127] + tmpQ2[62]*tmpFx[141] + tmpQ2[63]*tmpFx[155] + tmpQ2[64]*tmpFx[169];
tmpQ1[58] = + tmpQ2[52]*tmpFx[2] + tmpQ2[53]*tmpFx[16] + tmpQ2[54]*tmpFx[30] + tmpQ2[55]*tmpFx[44] + tmpQ2[56]*tmpFx[58] + tmpQ2[57]*tmpFx[72] + tmpQ2[58]*tmpFx[86] + tmpQ2[59]*tmpFx[100] + tmpQ2[60]*tmpFx[114] + tmpQ2[61]*tmpFx[128] + tmpQ2[62]*tmpFx[142] + tmpQ2[63]*tmpFx[156] + tmpQ2[64]*tmpFx[170];
tmpQ1[59] = + tmpQ2[52]*tmpFx[3] + tmpQ2[53]*tmpFx[17] + tmpQ2[54]*tmpFx[31] + tmpQ2[55]*tmpFx[45] + tmpQ2[56]*tmpFx[59] + tmpQ2[57]*tmpFx[73] + tmpQ2[58]*tmpFx[87] + tmpQ2[59]*tmpFx[101] + tmpQ2[60]*tmpFx[115] + tmpQ2[61]*tmpFx[129] + tmpQ2[62]*tmpFx[143] + tmpQ2[63]*tmpFx[157] + tmpQ2[64]*tmpFx[171];
tmpQ1[60] = + tmpQ2[52]*tmpFx[4] + tmpQ2[53]*tmpFx[18] + tmpQ2[54]*tmpFx[32] + tmpQ2[55]*tmpFx[46] + tmpQ2[56]*tmpFx[60] + tmpQ2[57]*tmpFx[74] + tmpQ2[58]*tmpFx[88] + tmpQ2[59]*tmpFx[102] + tmpQ2[60]*tmpFx[116] + tmpQ2[61]*tmpFx[130] + tmpQ2[62]*tmpFx[144] + tmpQ2[63]*tmpFx[158] + tmpQ2[64]*tmpFx[172];
tmpQ1[61] = + tmpQ2[52]*tmpFx[5] + tmpQ2[53]*tmpFx[19] + tmpQ2[54]*tmpFx[33] + tmpQ2[55]*tmpFx[47] + tmpQ2[56]*tmpFx[61] + tmpQ2[57]*tmpFx[75] + tmpQ2[58]*tmpFx[89] + tmpQ2[59]*tmpFx[103] + tmpQ2[60]*tmpFx[117] + tmpQ2[61]*tmpFx[131] + tmpQ2[62]*tmpFx[145] + tmpQ2[63]*tmpFx[159] + tmpQ2[64]*tmpFx[173];
tmpQ1[62] = + tmpQ2[52]*tmpFx[6] + tmpQ2[53]*tmpFx[20] + tmpQ2[54]*tmpFx[34] + tmpQ2[55]*tmpFx[48] + tmpQ2[56]*tmpFx[62] + tmpQ2[57]*tmpFx[76] + tmpQ2[58]*tmpFx[90] + tmpQ2[59]*tmpFx[104] + tmpQ2[60]*tmpFx[118] + tmpQ2[61]*tmpFx[132] + tmpQ2[62]*tmpFx[146] + tmpQ2[63]*tmpFx[160] + tmpQ2[64]*tmpFx[174];
tmpQ1[63] = + tmpQ2[52]*tmpFx[7] + tmpQ2[53]*tmpFx[21] + tmpQ2[54]*tmpFx[35] + tmpQ2[55]*tmpFx[49] + tmpQ2[56]*tmpFx[63] + tmpQ2[57]*tmpFx[77] + tmpQ2[58]*tmpFx[91] + tmpQ2[59]*tmpFx[105] + tmpQ2[60]*tmpFx[119] + tmpQ2[61]*tmpFx[133] + tmpQ2[62]*tmpFx[147] + tmpQ2[63]*tmpFx[161] + tmpQ2[64]*tmpFx[175];
tmpQ1[64] = + tmpQ2[52]*tmpFx[8] + tmpQ2[53]*tmpFx[22] + tmpQ2[54]*tmpFx[36] + tmpQ2[55]*tmpFx[50] + tmpQ2[56]*tmpFx[64] + tmpQ2[57]*tmpFx[78] + tmpQ2[58]*tmpFx[92] + tmpQ2[59]*tmpFx[106] + tmpQ2[60]*tmpFx[120] + tmpQ2[61]*tmpFx[134] + tmpQ2[62]*tmpFx[148] + tmpQ2[63]*tmpFx[162] + tmpQ2[64]*tmpFx[176];
tmpQ1[65] = + tmpQ2[52]*tmpFx[9] + tmpQ2[53]*tmpFx[23] + tmpQ2[54]*tmpFx[37] + tmpQ2[55]*tmpFx[51] + tmpQ2[56]*tmpFx[65] + tmpQ2[57]*tmpFx[79] + tmpQ2[58]*tmpFx[93] + tmpQ2[59]*tmpFx[107] + tmpQ2[60]*tmpFx[121] + tmpQ2[61]*tmpFx[135] + tmpQ2[62]*tmpFx[149] + tmpQ2[63]*tmpFx[163] + tmpQ2[64]*tmpFx[177];
tmpQ1[66] = + tmpQ2[52]*tmpFx[10] + tmpQ2[53]*tmpFx[24] + tmpQ2[54]*tmpFx[38] + tmpQ2[55]*tmpFx[52] + tmpQ2[56]*tmpFx[66] + tmpQ2[57]*tmpFx[80] + tmpQ2[58]*tmpFx[94] + tmpQ2[59]*tmpFx[108] + tmpQ2[60]*tmpFx[122] + tmpQ2[61]*tmpFx[136] + tmpQ2[62]*tmpFx[150] + tmpQ2[63]*tmpFx[164] + tmpQ2[64]*tmpFx[178];
tmpQ1[67] = + tmpQ2[52]*tmpFx[11] + tmpQ2[53]*tmpFx[25] + tmpQ2[54]*tmpFx[39] + tmpQ2[55]*tmpFx[53] + tmpQ2[56]*tmpFx[67] + tmpQ2[57]*tmpFx[81] + tmpQ2[58]*tmpFx[95] + tmpQ2[59]*tmpFx[109] + tmpQ2[60]*tmpFx[123] + tmpQ2[61]*tmpFx[137] + tmpQ2[62]*tmpFx[151] + tmpQ2[63]*tmpFx[165] + tmpQ2[64]*tmpFx[179];
tmpQ1[68] = + tmpQ2[52]*tmpFx[12] + tmpQ2[53]*tmpFx[26] + tmpQ2[54]*tmpFx[40] + tmpQ2[55]*tmpFx[54] + tmpQ2[56]*tmpFx[68] + tmpQ2[57]*tmpFx[82] + tmpQ2[58]*tmpFx[96] + tmpQ2[59]*tmpFx[110] + tmpQ2[60]*tmpFx[124] + tmpQ2[61]*tmpFx[138] + tmpQ2[62]*tmpFx[152] + tmpQ2[63]*tmpFx[166] + tmpQ2[64]*tmpFx[180];
tmpQ1[69] = + tmpQ2[52]*tmpFx[13] + tmpQ2[53]*tmpFx[27] + tmpQ2[54]*tmpFx[41] + tmpQ2[55]*tmpFx[55] + tmpQ2[56]*tmpFx[69] + tmpQ2[57]*tmpFx[83] + tmpQ2[58]*tmpFx[97] + tmpQ2[59]*tmpFx[111] + tmpQ2[60]*tmpFx[125] + tmpQ2[61]*tmpFx[139] + tmpQ2[62]*tmpFx[153] + tmpQ2[63]*tmpFx[167] + tmpQ2[64]*tmpFx[181];
tmpQ1[70] = + tmpQ2[65]*tmpFx[0] + tmpQ2[66]*tmpFx[14] + tmpQ2[67]*tmpFx[28] + tmpQ2[68]*tmpFx[42] + tmpQ2[69]*tmpFx[56] + tmpQ2[70]*tmpFx[70] + tmpQ2[71]*tmpFx[84] + tmpQ2[72]*tmpFx[98] + tmpQ2[73]*tmpFx[112] + tmpQ2[74]*tmpFx[126] + tmpQ2[75]*tmpFx[140] + tmpQ2[76]*tmpFx[154] + tmpQ2[77]*tmpFx[168];
tmpQ1[71] = + tmpQ2[65]*tmpFx[1] + tmpQ2[66]*tmpFx[15] + tmpQ2[67]*tmpFx[29] + tmpQ2[68]*tmpFx[43] + tmpQ2[69]*tmpFx[57] + tmpQ2[70]*tmpFx[71] + tmpQ2[71]*tmpFx[85] + tmpQ2[72]*tmpFx[99] + tmpQ2[73]*tmpFx[113] + tmpQ2[74]*tmpFx[127] + tmpQ2[75]*tmpFx[141] + tmpQ2[76]*tmpFx[155] + tmpQ2[77]*tmpFx[169];
tmpQ1[72] = + tmpQ2[65]*tmpFx[2] + tmpQ2[66]*tmpFx[16] + tmpQ2[67]*tmpFx[30] + tmpQ2[68]*tmpFx[44] + tmpQ2[69]*tmpFx[58] + tmpQ2[70]*tmpFx[72] + tmpQ2[71]*tmpFx[86] + tmpQ2[72]*tmpFx[100] + tmpQ2[73]*tmpFx[114] + tmpQ2[74]*tmpFx[128] + tmpQ2[75]*tmpFx[142] + tmpQ2[76]*tmpFx[156] + tmpQ2[77]*tmpFx[170];
tmpQ1[73] = + tmpQ2[65]*tmpFx[3] + tmpQ2[66]*tmpFx[17] + tmpQ2[67]*tmpFx[31] + tmpQ2[68]*tmpFx[45] + tmpQ2[69]*tmpFx[59] + tmpQ2[70]*tmpFx[73] + tmpQ2[71]*tmpFx[87] + tmpQ2[72]*tmpFx[101] + tmpQ2[73]*tmpFx[115] + tmpQ2[74]*tmpFx[129] + tmpQ2[75]*tmpFx[143] + tmpQ2[76]*tmpFx[157] + tmpQ2[77]*tmpFx[171];
tmpQ1[74] = + tmpQ2[65]*tmpFx[4] + tmpQ2[66]*tmpFx[18] + tmpQ2[67]*tmpFx[32] + tmpQ2[68]*tmpFx[46] + tmpQ2[69]*tmpFx[60] + tmpQ2[70]*tmpFx[74] + tmpQ2[71]*tmpFx[88] + tmpQ2[72]*tmpFx[102] + tmpQ2[73]*tmpFx[116] + tmpQ2[74]*tmpFx[130] + tmpQ2[75]*tmpFx[144] + tmpQ2[76]*tmpFx[158] + tmpQ2[77]*tmpFx[172];
tmpQ1[75] = + tmpQ2[65]*tmpFx[5] + tmpQ2[66]*tmpFx[19] + tmpQ2[67]*tmpFx[33] + tmpQ2[68]*tmpFx[47] + tmpQ2[69]*tmpFx[61] + tmpQ2[70]*tmpFx[75] + tmpQ2[71]*tmpFx[89] + tmpQ2[72]*tmpFx[103] + tmpQ2[73]*tmpFx[117] + tmpQ2[74]*tmpFx[131] + tmpQ2[75]*tmpFx[145] + tmpQ2[76]*tmpFx[159] + tmpQ2[77]*tmpFx[173];
tmpQ1[76] = + tmpQ2[65]*tmpFx[6] + tmpQ2[66]*tmpFx[20] + tmpQ2[67]*tmpFx[34] + tmpQ2[68]*tmpFx[48] + tmpQ2[69]*tmpFx[62] + tmpQ2[70]*tmpFx[76] + tmpQ2[71]*tmpFx[90] + tmpQ2[72]*tmpFx[104] + tmpQ2[73]*tmpFx[118] + tmpQ2[74]*tmpFx[132] + tmpQ2[75]*tmpFx[146] + tmpQ2[76]*tmpFx[160] + tmpQ2[77]*tmpFx[174];
tmpQ1[77] = + tmpQ2[65]*tmpFx[7] + tmpQ2[66]*tmpFx[21] + tmpQ2[67]*tmpFx[35] + tmpQ2[68]*tmpFx[49] + tmpQ2[69]*tmpFx[63] + tmpQ2[70]*tmpFx[77] + tmpQ2[71]*tmpFx[91] + tmpQ2[72]*tmpFx[105] + tmpQ2[73]*tmpFx[119] + tmpQ2[74]*tmpFx[133] + tmpQ2[75]*tmpFx[147] + tmpQ2[76]*tmpFx[161] + tmpQ2[77]*tmpFx[175];
tmpQ1[78] = + tmpQ2[65]*tmpFx[8] + tmpQ2[66]*tmpFx[22] + tmpQ2[67]*tmpFx[36] + tmpQ2[68]*tmpFx[50] + tmpQ2[69]*tmpFx[64] + tmpQ2[70]*tmpFx[78] + tmpQ2[71]*tmpFx[92] + tmpQ2[72]*tmpFx[106] + tmpQ2[73]*tmpFx[120] + tmpQ2[74]*tmpFx[134] + tmpQ2[75]*tmpFx[148] + tmpQ2[76]*tmpFx[162] + tmpQ2[77]*tmpFx[176];
tmpQ1[79] = + tmpQ2[65]*tmpFx[9] + tmpQ2[66]*tmpFx[23] + tmpQ2[67]*tmpFx[37] + tmpQ2[68]*tmpFx[51] + tmpQ2[69]*tmpFx[65] + tmpQ2[70]*tmpFx[79] + tmpQ2[71]*tmpFx[93] + tmpQ2[72]*tmpFx[107] + tmpQ2[73]*tmpFx[121] + tmpQ2[74]*tmpFx[135] + tmpQ2[75]*tmpFx[149] + tmpQ2[76]*tmpFx[163] + tmpQ2[77]*tmpFx[177];
tmpQ1[80] = + tmpQ2[65]*tmpFx[10] + tmpQ2[66]*tmpFx[24] + tmpQ2[67]*tmpFx[38] + tmpQ2[68]*tmpFx[52] + tmpQ2[69]*tmpFx[66] + tmpQ2[70]*tmpFx[80] + tmpQ2[71]*tmpFx[94] + tmpQ2[72]*tmpFx[108] + tmpQ2[73]*tmpFx[122] + tmpQ2[74]*tmpFx[136] + tmpQ2[75]*tmpFx[150] + tmpQ2[76]*tmpFx[164] + tmpQ2[77]*tmpFx[178];
tmpQ1[81] = + tmpQ2[65]*tmpFx[11] + tmpQ2[66]*tmpFx[25] + tmpQ2[67]*tmpFx[39] + tmpQ2[68]*tmpFx[53] + tmpQ2[69]*tmpFx[67] + tmpQ2[70]*tmpFx[81] + tmpQ2[71]*tmpFx[95] + tmpQ2[72]*tmpFx[109] + tmpQ2[73]*tmpFx[123] + tmpQ2[74]*tmpFx[137] + tmpQ2[75]*tmpFx[151] + tmpQ2[76]*tmpFx[165] + tmpQ2[77]*tmpFx[179];
tmpQ1[82] = + tmpQ2[65]*tmpFx[12] + tmpQ2[66]*tmpFx[26] + tmpQ2[67]*tmpFx[40] + tmpQ2[68]*tmpFx[54] + tmpQ2[69]*tmpFx[68] + tmpQ2[70]*tmpFx[82] + tmpQ2[71]*tmpFx[96] + tmpQ2[72]*tmpFx[110] + tmpQ2[73]*tmpFx[124] + tmpQ2[74]*tmpFx[138] + tmpQ2[75]*tmpFx[152] + tmpQ2[76]*tmpFx[166] + tmpQ2[77]*tmpFx[180];
tmpQ1[83] = + tmpQ2[65]*tmpFx[13] + tmpQ2[66]*tmpFx[27] + tmpQ2[67]*tmpFx[41] + tmpQ2[68]*tmpFx[55] + tmpQ2[69]*tmpFx[69] + tmpQ2[70]*tmpFx[83] + tmpQ2[71]*tmpFx[97] + tmpQ2[72]*tmpFx[111] + tmpQ2[73]*tmpFx[125] + tmpQ2[74]*tmpFx[139] + tmpQ2[75]*tmpFx[153] + tmpQ2[76]*tmpFx[167] + tmpQ2[77]*tmpFx[181];
tmpQ1[84] = + tmpQ2[78]*tmpFx[0] + tmpQ2[79]*tmpFx[14] + tmpQ2[80]*tmpFx[28] + tmpQ2[81]*tmpFx[42] + tmpQ2[82]*tmpFx[56] + tmpQ2[83]*tmpFx[70] + tmpQ2[84]*tmpFx[84] + tmpQ2[85]*tmpFx[98] + tmpQ2[86]*tmpFx[112] + tmpQ2[87]*tmpFx[126] + tmpQ2[88]*tmpFx[140] + tmpQ2[89]*tmpFx[154] + tmpQ2[90]*tmpFx[168];
tmpQ1[85] = + tmpQ2[78]*tmpFx[1] + tmpQ2[79]*tmpFx[15] + tmpQ2[80]*tmpFx[29] + tmpQ2[81]*tmpFx[43] + tmpQ2[82]*tmpFx[57] + tmpQ2[83]*tmpFx[71] + tmpQ2[84]*tmpFx[85] + tmpQ2[85]*tmpFx[99] + tmpQ2[86]*tmpFx[113] + tmpQ2[87]*tmpFx[127] + tmpQ2[88]*tmpFx[141] + tmpQ2[89]*tmpFx[155] + tmpQ2[90]*tmpFx[169];
tmpQ1[86] = + tmpQ2[78]*tmpFx[2] + tmpQ2[79]*tmpFx[16] + tmpQ2[80]*tmpFx[30] + tmpQ2[81]*tmpFx[44] + tmpQ2[82]*tmpFx[58] + tmpQ2[83]*tmpFx[72] + tmpQ2[84]*tmpFx[86] + tmpQ2[85]*tmpFx[100] + tmpQ2[86]*tmpFx[114] + tmpQ2[87]*tmpFx[128] + tmpQ2[88]*tmpFx[142] + tmpQ2[89]*tmpFx[156] + tmpQ2[90]*tmpFx[170];
tmpQ1[87] = + tmpQ2[78]*tmpFx[3] + tmpQ2[79]*tmpFx[17] + tmpQ2[80]*tmpFx[31] + tmpQ2[81]*tmpFx[45] + tmpQ2[82]*tmpFx[59] + tmpQ2[83]*tmpFx[73] + tmpQ2[84]*tmpFx[87] + tmpQ2[85]*tmpFx[101] + tmpQ2[86]*tmpFx[115] + tmpQ2[87]*tmpFx[129] + tmpQ2[88]*tmpFx[143] + tmpQ2[89]*tmpFx[157] + tmpQ2[90]*tmpFx[171];
tmpQ1[88] = + tmpQ2[78]*tmpFx[4] + tmpQ2[79]*tmpFx[18] + tmpQ2[80]*tmpFx[32] + tmpQ2[81]*tmpFx[46] + tmpQ2[82]*tmpFx[60] + tmpQ2[83]*tmpFx[74] + tmpQ2[84]*tmpFx[88] + tmpQ2[85]*tmpFx[102] + tmpQ2[86]*tmpFx[116] + tmpQ2[87]*tmpFx[130] + tmpQ2[88]*tmpFx[144] + tmpQ2[89]*tmpFx[158] + tmpQ2[90]*tmpFx[172];
tmpQ1[89] = + tmpQ2[78]*tmpFx[5] + tmpQ2[79]*tmpFx[19] + tmpQ2[80]*tmpFx[33] + tmpQ2[81]*tmpFx[47] + tmpQ2[82]*tmpFx[61] + tmpQ2[83]*tmpFx[75] + tmpQ2[84]*tmpFx[89] + tmpQ2[85]*tmpFx[103] + tmpQ2[86]*tmpFx[117] + tmpQ2[87]*tmpFx[131] + tmpQ2[88]*tmpFx[145] + tmpQ2[89]*tmpFx[159] + tmpQ2[90]*tmpFx[173];
tmpQ1[90] = + tmpQ2[78]*tmpFx[6] + tmpQ2[79]*tmpFx[20] + tmpQ2[80]*tmpFx[34] + tmpQ2[81]*tmpFx[48] + tmpQ2[82]*tmpFx[62] + tmpQ2[83]*tmpFx[76] + tmpQ2[84]*tmpFx[90] + tmpQ2[85]*tmpFx[104] + tmpQ2[86]*tmpFx[118] + tmpQ2[87]*tmpFx[132] + tmpQ2[88]*tmpFx[146] + tmpQ2[89]*tmpFx[160] + tmpQ2[90]*tmpFx[174];
tmpQ1[91] = + tmpQ2[78]*tmpFx[7] + tmpQ2[79]*tmpFx[21] + tmpQ2[80]*tmpFx[35] + tmpQ2[81]*tmpFx[49] + tmpQ2[82]*tmpFx[63] + tmpQ2[83]*tmpFx[77] + tmpQ2[84]*tmpFx[91] + tmpQ2[85]*tmpFx[105] + tmpQ2[86]*tmpFx[119] + tmpQ2[87]*tmpFx[133] + tmpQ2[88]*tmpFx[147] + tmpQ2[89]*tmpFx[161] + tmpQ2[90]*tmpFx[175];
tmpQ1[92] = + tmpQ2[78]*tmpFx[8] + tmpQ2[79]*tmpFx[22] + tmpQ2[80]*tmpFx[36] + tmpQ2[81]*tmpFx[50] + tmpQ2[82]*tmpFx[64] + tmpQ2[83]*tmpFx[78] + tmpQ2[84]*tmpFx[92] + tmpQ2[85]*tmpFx[106] + tmpQ2[86]*tmpFx[120] + tmpQ2[87]*tmpFx[134] + tmpQ2[88]*tmpFx[148] + tmpQ2[89]*tmpFx[162] + tmpQ2[90]*tmpFx[176];
tmpQ1[93] = + tmpQ2[78]*tmpFx[9] + tmpQ2[79]*tmpFx[23] + tmpQ2[80]*tmpFx[37] + tmpQ2[81]*tmpFx[51] + tmpQ2[82]*tmpFx[65] + tmpQ2[83]*tmpFx[79] + tmpQ2[84]*tmpFx[93] + tmpQ2[85]*tmpFx[107] + tmpQ2[86]*tmpFx[121] + tmpQ2[87]*tmpFx[135] + tmpQ2[88]*tmpFx[149] + tmpQ2[89]*tmpFx[163] + tmpQ2[90]*tmpFx[177];
tmpQ1[94] = + tmpQ2[78]*tmpFx[10] + tmpQ2[79]*tmpFx[24] + tmpQ2[80]*tmpFx[38] + tmpQ2[81]*tmpFx[52] + tmpQ2[82]*tmpFx[66] + tmpQ2[83]*tmpFx[80] + tmpQ2[84]*tmpFx[94] + tmpQ2[85]*tmpFx[108] + tmpQ2[86]*tmpFx[122] + tmpQ2[87]*tmpFx[136] + tmpQ2[88]*tmpFx[150] + tmpQ2[89]*tmpFx[164] + tmpQ2[90]*tmpFx[178];
tmpQ1[95] = + tmpQ2[78]*tmpFx[11] + tmpQ2[79]*tmpFx[25] + tmpQ2[80]*tmpFx[39] + tmpQ2[81]*tmpFx[53] + tmpQ2[82]*tmpFx[67] + tmpQ2[83]*tmpFx[81] + tmpQ2[84]*tmpFx[95] + tmpQ2[85]*tmpFx[109] + tmpQ2[86]*tmpFx[123] + tmpQ2[87]*tmpFx[137] + tmpQ2[88]*tmpFx[151] + tmpQ2[89]*tmpFx[165] + tmpQ2[90]*tmpFx[179];
tmpQ1[96] = + tmpQ2[78]*tmpFx[12] + tmpQ2[79]*tmpFx[26] + tmpQ2[80]*tmpFx[40] + tmpQ2[81]*tmpFx[54] + tmpQ2[82]*tmpFx[68] + tmpQ2[83]*tmpFx[82] + tmpQ2[84]*tmpFx[96] + tmpQ2[85]*tmpFx[110] + tmpQ2[86]*tmpFx[124] + tmpQ2[87]*tmpFx[138] + tmpQ2[88]*tmpFx[152] + tmpQ2[89]*tmpFx[166] + tmpQ2[90]*tmpFx[180];
tmpQ1[97] = + tmpQ2[78]*tmpFx[13] + tmpQ2[79]*tmpFx[27] + tmpQ2[80]*tmpFx[41] + tmpQ2[81]*tmpFx[55] + tmpQ2[82]*tmpFx[69] + tmpQ2[83]*tmpFx[83] + tmpQ2[84]*tmpFx[97] + tmpQ2[85]*tmpFx[111] + tmpQ2[86]*tmpFx[125] + tmpQ2[87]*tmpFx[139] + tmpQ2[88]*tmpFx[153] + tmpQ2[89]*tmpFx[167] + tmpQ2[90]*tmpFx[181];
tmpQ1[98] = + tmpQ2[91]*tmpFx[0] + tmpQ2[92]*tmpFx[14] + tmpQ2[93]*tmpFx[28] + tmpQ2[94]*tmpFx[42] + tmpQ2[95]*tmpFx[56] + tmpQ2[96]*tmpFx[70] + tmpQ2[97]*tmpFx[84] + tmpQ2[98]*tmpFx[98] + tmpQ2[99]*tmpFx[112] + tmpQ2[100]*tmpFx[126] + tmpQ2[101]*tmpFx[140] + tmpQ2[102]*tmpFx[154] + tmpQ2[103]*tmpFx[168];
tmpQ1[99] = + tmpQ2[91]*tmpFx[1] + tmpQ2[92]*tmpFx[15] + tmpQ2[93]*tmpFx[29] + tmpQ2[94]*tmpFx[43] + tmpQ2[95]*tmpFx[57] + tmpQ2[96]*tmpFx[71] + tmpQ2[97]*tmpFx[85] + tmpQ2[98]*tmpFx[99] + tmpQ2[99]*tmpFx[113] + tmpQ2[100]*tmpFx[127] + tmpQ2[101]*tmpFx[141] + tmpQ2[102]*tmpFx[155] + tmpQ2[103]*tmpFx[169];
tmpQ1[100] = + tmpQ2[91]*tmpFx[2] + tmpQ2[92]*tmpFx[16] + tmpQ2[93]*tmpFx[30] + tmpQ2[94]*tmpFx[44] + tmpQ2[95]*tmpFx[58] + tmpQ2[96]*tmpFx[72] + tmpQ2[97]*tmpFx[86] + tmpQ2[98]*tmpFx[100] + tmpQ2[99]*tmpFx[114] + tmpQ2[100]*tmpFx[128] + tmpQ2[101]*tmpFx[142] + tmpQ2[102]*tmpFx[156] + tmpQ2[103]*tmpFx[170];
tmpQ1[101] = + tmpQ2[91]*tmpFx[3] + tmpQ2[92]*tmpFx[17] + tmpQ2[93]*tmpFx[31] + tmpQ2[94]*tmpFx[45] + tmpQ2[95]*tmpFx[59] + tmpQ2[96]*tmpFx[73] + tmpQ2[97]*tmpFx[87] + tmpQ2[98]*tmpFx[101] + tmpQ2[99]*tmpFx[115] + tmpQ2[100]*tmpFx[129] + tmpQ2[101]*tmpFx[143] + tmpQ2[102]*tmpFx[157] + tmpQ2[103]*tmpFx[171];
tmpQ1[102] = + tmpQ2[91]*tmpFx[4] + tmpQ2[92]*tmpFx[18] + tmpQ2[93]*tmpFx[32] + tmpQ2[94]*tmpFx[46] + tmpQ2[95]*tmpFx[60] + tmpQ2[96]*tmpFx[74] + tmpQ2[97]*tmpFx[88] + tmpQ2[98]*tmpFx[102] + tmpQ2[99]*tmpFx[116] + tmpQ2[100]*tmpFx[130] + tmpQ2[101]*tmpFx[144] + tmpQ2[102]*tmpFx[158] + tmpQ2[103]*tmpFx[172];
tmpQ1[103] = + tmpQ2[91]*tmpFx[5] + tmpQ2[92]*tmpFx[19] + tmpQ2[93]*tmpFx[33] + tmpQ2[94]*tmpFx[47] + tmpQ2[95]*tmpFx[61] + tmpQ2[96]*tmpFx[75] + tmpQ2[97]*tmpFx[89] + tmpQ2[98]*tmpFx[103] + tmpQ2[99]*tmpFx[117] + tmpQ2[100]*tmpFx[131] + tmpQ2[101]*tmpFx[145] + tmpQ2[102]*tmpFx[159] + tmpQ2[103]*tmpFx[173];
tmpQ1[104] = + tmpQ2[91]*tmpFx[6] + tmpQ2[92]*tmpFx[20] + tmpQ2[93]*tmpFx[34] + tmpQ2[94]*tmpFx[48] + tmpQ2[95]*tmpFx[62] + tmpQ2[96]*tmpFx[76] + tmpQ2[97]*tmpFx[90] + tmpQ2[98]*tmpFx[104] + tmpQ2[99]*tmpFx[118] + tmpQ2[100]*tmpFx[132] + tmpQ2[101]*tmpFx[146] + tmpQ2[102]*tmpFx[160] + tmpQ2[103]*tmpFx[174];
tmpQ1[105] = + tmpQ2[91]*tmpFx[7] + tmpQ2[92]*tmpFx[21] + tmpQ2[93]*tmpFx[35] + tmpQ2[94]*tmpFx[49] + tmpQ2[95]*tmpFx[63] + tmpQ2[96]*tmpFx[77] + tmpQ2[97]*tmpFx[91] + tmpQ2[98]*tmpFx[105] + tmpQ2[99]*tmpFx[119] + tmpQ2[100]*tmpFx[133] + tmpQ2[101]*tmpFx[147] + tmpQ2[102]*tmpFx[161] + tmpQ2[103]*tmpFx[175];
tmpQ1[106] = + tmpQ2[91]*tmpFx[8] + tmpQ2[92]*tmpFx[22] + tmpQ2[93]*tmpFx[36] + tmpQ2[94]*tmpFx[50] + tmpQ2[95]*tmpFx[64] + tmpQ2[96]*tmpFx[78] + tmpQ2[97]*tmpFx[92] + tmpQ2[98]*tmpFx[106] + tmpQ2[99]*tmpFx[120] + tmpQ2[100]*tmpFx[134] + tmpQ2[101]*tmpFx[148] + tmpQ2[102]*tmpFx[162] + tmpQ2[103]*tmpFx[176];
tmpQ1[107] = + tmpQ2[91]*tmpFx[9] + tmpQ2[92]*tmpFx[23] + tmpQ2[93]*tmpFx[37] + tmpQ2[94]*tmpFx[51] + tmpQ2[95]*tmpFx[65] + tmpQ2[96]*tmpFx[79] + tmpQ2[97]*tmpFx[93] + tmpQ2[98]*tmpFx[107] + tmpQ2[99]*tmpFx[121] + tmpQ2[100]*tmpFx[135] + tmpQ2[101]*tmpFx[149] + tmpQ2[102]*tmpFx[163] + tmpQ2[103]*tmpFx[177];
tmpQ1[108] = + tmpQ2[91]*tmpFx[10] + tmpQ2[92]*tmpFx[24] + tmpQ2[93]*tmpFx[38] + tmpQ2[94]*tmpFx[52] + tmpQ2[95]*tmpFx[66] + tmpQ2[96]*tmpFx[80] + tmpQ2[97]*tmpFx[94] + tmpQ2[98]*tmpFx[108] + tmpQ2[99]*tmpFx[122] + tmpQ2[100]*tmpFx[136] + tmpQ2[101]*tmpFx[150] + tmpQ2[102]*tmpFx[164] + tmpQ2[103]*tmpFx[178];
tmpQ1[109] = + tmpQ2[91]*tmpFx[11] + tmpQ2[92]*tmpFx[25] + tmpQ2[93]*tmpFx[39] + tmpQ2[94]*tmpFx[53] + tmpQ2[95]*tmpFx[67] + tmpQ2[96]*tmpFx[81] + tmpQ2[97]*tmpFx[95] + tmpQ2[98]*tmpFx[109] + tmpQ2[99]*tmpFx[123] + tmpQ2[100]*tmpFx[137] + tmpQ2[101]*tmpFx[151] + tmpQ2[102]*tmpFx[165] + tmpQ2[103]*tmpFx[179];
tmpQ1[110] = + tmpQ2[91]*tmpFx[12] + tmpQ2[92]*tmpFx[26] + tmpQ2[93]*tmpFx[40] + tmpQ2[94]*tmpFx[54] + tmpQ2[95]*tmpFx[68] + tmpQ2[96]*tmpFx[82] + tmpQ2[97]*tmpFx[96] + tmpQ2[98]*tmpFx[110] + tmpQ2[99]*tmpFx[124] + tmpQ2[100]*tmpFx[138] + tmpQ2[101]*tmpFx[152] + tmpQ2[102]*tmpFx[166] + tmpQ2[103]*tmpFx[180];
tmpQ1[111] = + tmpQ2[91]*tmpFx[13] + tmpQ2[92]*tmpFx[27] + tmpQ2[93]*tmpFx[41] + tmpQ2[94]*tmpFx[55] + tmpQ2[95]*tmpFx[69] + tmpQ2[96]*tmpFx[83] + tmpQ2[97]*tmpFx[97] + tmpQ2[98]*tmpFx[111] + tmpQ2[99]*tmpFx[125] + tmpQ2[100]*tmpFx[139] + tmpQ2[101]*tmpFx[153] + tmpQ2[102]*tmpFx[167] + tmpQ2[103]*tmpFx[181];
tmpQ1[112] = + tmpQ2[104]*tmpFx[0] + tmpQ2[105]*tmpFx[14] + tmpQ2[106]*tmpFx[28] + tmpQ2[107]*tmpFx[42] + tmpQ2[108]*tmpFx[56] + tmpQ2[109]*tmpFx[70] + tmpQ2[110]*tmpFx[84] + tmpQ2[111]*tmpFx[98] + tmpQ2[112]*tmpFx[112] + tmpQ2[113]*tmpFx[126] + tmpQ2[114]*tmpFx[140] + tmpQ2[115]*tmpFx[154] + tmpQ2[116]*tmpFx[168];
tmpQ1[113] = + tmpQ2[104]*tmpFx[1] + tmpQ2[105]*tmpFx[15] + tmpQ2[106]*tmpFx[29] + tmpQ2[107]*tmpFx[43] + tmpQ2[108]*tmpFx[57] + tmpQ2[109]*tmpFx[71] + tmpQ2[110]*tmpFx[85] + tmpQ2[111]*tmpFx[99] + tmpQ2[112]*tmpFx[113] + tmpQ2[113]*tmpFx[127] + tmpQ2[114]*tmpFx[141] + tmpQ2[115]*tmpFx[155] + tmpQ2[116]*tmpFx[169];
tmpQ1[114] = + tmpQ2[104]*tmpFx[2] + tmpQ2[105]*tmpFx[16] + tmpQ2[106]*tmpFx[30] + tmpQ2[107]*tmpFx[44] + tmpQ2[108]*tmpFx[58] + tmpQ2[109]*tmpFx[72] + tmpQ2[110]*tmpFx[86] + tmpQ2[111]*tmpFx[100] + tmpQ2[112]*tmpFx[114] + tmpQ2[113]*tmpFx[128] + tmpQ2[114]*tmpFx[142] + tmpQ2[115]*tmpFx[156] + tmpQ2[116]*tmpFx[170];
tmpQ1[115] = + tmpQ2[104]*tmpFx[3] + tmpQ2[105]*tmpFx[17] + tmpQ2[106]*tmpFx[31] + tmpQ2[107]*tmpFx[45] + tmpQ2[108]*tmpFx[59] + tmpQ2[109]*tmpFx[73] + tmpQ2[110]*tmpFx[87] + tmpQ2[111]*tmpFx[101] + tmpQ2[112]*tmpFx[115] + tmpQ2[113]*tmpFx[129] + tmpQ2[114]*tmpFx[143] + tmpQ2[115]*tmpFx[157] + tmpQ2[116]*tmpFx[171];
tmpQ1[116] = + tmpQ2[104]*tmpFx[4] + tmpQ2[105]*tmpFx[18] + tmpQ2[106]*tmpFx[32] + tmpQ2[107]*tmpFx[46] + tmpQ2[108]*tmpFx[60] + tmpQ2[109]*tmpFx[74] + tmpQ2[110]*tmpFx[88] + tmpQ2[111]*tmpFx[102] + tmpQ2[112]*tmpFx[116] + tmpQ2[113]*tmpFx[130] + tmpQ2[114]*tmpFx[144] + tmpQ2[115]*tmpFx[158] + tmpQ2[116]*tmpFx[172];
tmpQ1[117] = + tmpQ2[104]*tmpFx[5] + tmpQ2[105]*tmpFx[19] + tmpQ2[106]*tmpFx[33] + tmpQ2[107]*tmpFx[47] + tmpQ2[108]*tmpFx[61] + tmpQ2[109]*tmpFx[75] + tmpQ2[110]*tmpFx[89] + tmpQ2[111]*tmpFx[103] + tmpQ2[112]*tmpFx[117] + tmpQ2[113]*tmpFx[131] + tmpQ2[114]*tmpFx[145] + tmpQ2[115]*tmpFx[159] + tmpQ2[116]*tmpFx[173];
tmpQ1[118] = + tmpQ2[104]*tmpFx[6] + tmpQ2[105]*tmpFx[20] + tmpQ2[106]*tmpFx[34] + tmpQ2[107]*tmpFx[48] + tmpQ2[108]*tmpFx[62] + tmpQ2[109]*tmpFx[76] + tmpQ2[110]*tmpFx[90] + tmpQ2[111]*tmpFx[104] + tmpQ2[112]*tmpFx[118] + tmpQ2[113]*tmpFx[132] + tmpQ2[114]*tmpFx[146] + tmpQ2[115]*tmpFx[160] + tmpQ2[116]*tmpFx[174];
tmpQ1[119] = + tmpQ2[104]*tmpFx[7] + tmpQ2[105]*tmpFx[21] + tmpQ2[106]*tmpFx[35] + tmpQ2[107]*tmpFx[49] + tmpQ2[108]*tmpFx[63] + tmpQ2[109]*tmpFx[77] + tmpQ2[110]*tmpFx[91] + tmpQ2[111]*tmpFx[105] + tmpQ2[112]*tmpFx[119] + tmpQ2[113]*tmpFx[133] + tmpQ2[114]*tmpFx[147] + tmpQ2[115]*tmpFx[161] + tmpQ2[116]*tmpFx[175];
tmpQ1[120] = + tmpQ2[104]*tmpFx[8] + tmpQ2[105]*tmpFx[22] + tmpQ2[106]*tmpFx[36] + tmpQ2[107]*tmpFx[50] + tmpQ2[108]*tmpFx[64] + tmpQ2[109]*tmpFx[78] + tmpQ2[110]*tmpFx[92] + tmpQ2[111]*tmpFx[106] + tmpQ2[112]*tmpFx[120] + tmpQ2[113]*tmpFx[134] + tmpQ2[114]*tmpFx[148] + tmpQ2[115]*tmpFx[162] + tmpQ2[116]*tmpFx[176];
tmpQ1[121] = + tmpQ2[104]*tmpFx[9] + tmpQ2[105]*tmpFx[23] + tmpQ2[106]*tmpFx[37] + tmpQ2[107]*tmpFx[51] + tmpQ2[108]*tmpFx[65] + tmpQ2[109]*tmpFx[79] + tmpQ2[110]*tmpFx[93] + tmpQ2[111]*tmpFx[107] + tmpQ2[112]*tmpFx[121] + tmpQ2[113]*tmpFx[135] + tmpQ2[114]*tmpFx[149] + tmpQ2[115]*tmpFx[163] + tmpQ2[116]*tmpFx[177];
tmpQ1[122] = + tmpQ2[104]*tmpFx[10] + tmpQ2[105]*tmpFx[24] + tmpQ2[106]*tmpFx[38] + tmpQ2[107]*tmpFx[52] + tmpQ2[108]*tmpFx[66] + tmpQ2[109]*tmpFx[80] + tmpQ2[110]*tmpFx[94] + tmpQ2[111]*tmpFx[108] + tmpQ2[112]*tmpFx[122] + tmpQ2[113]*tmpFx[136] + tmpQ2[114]*tmpFx[150] + tmpQ2[115]*tmpFx[164] + tmpQ2[116]*tmpFx[178];
tmpQ1[123] = + tmpQ2[104]*tmpFx[11] + tmpQ2[105]*tmpFx[25] + tmpQ2[106]*tmpFx[39] + tmpQ2[107]*tmpFx[53] + tmpQ2[108]*tmpFx[67] + tmpQ2[109]*tmpFx[81] + tmpQ2[110]*tmpFx[95] + tmpQ2[111]*tmpFx[109] + tmpQ2[112]*tmpFx[123] + tmpQ2[113]*tmpFx[137] + tmpQ2[114]*tmpFx[151] + tmpQ2[115]*tmpFx[165] + tmpQ2[116]*tmpFx[179];
tmpQ1[124] = + tmpQ2[104]*tmpFx[12] + tmpQ2[105]*tmpFx[26] + tmpQ2[106]*tmpFx[40] + tmpQ2[107]*tmpFx[54] + tmpQ2[108]*tmpFx[68] + tmpQ2[109]*tmpFx[82] + tmpQ2[110]*tmpFx[96] + tmpQ2[111]*tmpFx[110] + tmpQ2[112]*tmpFx[124] + tmpQ2[113]*tmpFx[138] + tmpQ2[114]*tmpFx[152] + tmpQ2[115]*tmpFx[166] + tmpQ2[116]*tmpFx[180];
tmpQ1[125] = + tmpQ2[104]*tmpFx[13] + tmpQ2[105]*tmpFx[27] + tmpQ2[106]*tmpFx[41] + tmpQ2[107]*tmpFx[55] + tmpQ2[108]*tmpFx[69] + tmpQ2[109]*tmpFx[83] + tmpQ2[110]*tmpFx[97] + tmpQ2[111]*tmpFx[111] + tmpQ2[112]*tmpFx[125] + tmpQ2[113]*tmpFx[139] + tmpQ2[114]*tmpFx[153] + tmpQ2[115]*tmpFx[167] + tmpQ2[116]*tmpFx[181];
tmpQ1[126] = + tmpQ2[117]*tmpFx[0] + tmpQ2[118]*tmpFx[14] + tmpQ2[119]*tmpFx[28] + tmpQ2[120]*tmpFx[42] + tmpQ2[121]*tmpFx[56] + tmpQ2[122]*tmpFx[70] + tmpQ2[123]*tmpFx[84] + tmpQ2[124]*tmpFx[98] + tmpQ2[125]*tmpFx[112] + tmpQ2[126]*tmpFx[126] + tmpQ2[127]*tmpFx[140] + tmpQ2[128]*tmpFx[154] + tmpQ2[129]*tmpFx[168];
tmpQ1[127] = + tmpQ2[117]*tmpFx[1] + tmpQ2[118]*tmpFx[15] + tmpQ2[119]*tmpFx[29] + tmpQ2[120]*tmpFx[43] + tmpQ2[121]*tmpFx[57] + tmpQ2[122]*tmpFx[71] + tmpQ2[123]*tmpFx[85] + tmpQ2[124]*tmpFx[99] + tmpQ2[125]*tmpFx[113] + tmpQ2[126]*tmpFx[127] + tmpQ2[127]*tmpFx[141] + tmpQ2[128]*tmpFx[155] + tmpQ2[129]*tmpFx[169];
tmpQ1[128] = + tmpQ2[117]*tmpFx[2] + tmpQ2[118]*tmpFx[16] + tmpQ2[119]*tmpFx[30] + tmpQ2[120]*tmpFx[44] + tmpQ2[121]*tmpFx[58] + tmpQ2[122]*tmpFx[72] + tmpQ2[123]*tmpFx[86] + tmpQ2[124]*tmpFx[100] + tmpQ2[125]*tmpFx[114] + tmpQ2[126]*tmpFx[128] + tmpQ2[127]*tmpFx[142] + tmpQ2[128]*tmpFx[156] + tmpQ2[129]*tmpFx[170];
tmpQ1[129] = + tmpQ2[117]*tmpFx[3] + tmpQ2[118]*tmpFx[17] + tmpQ2[119]*tmpFx[31] + tmpQ2[120]*tmpFx[45] + tmpQ2[121]*tmpFx[59] + tmpQ2[122]*tmpFx[73] + tmpQ2[123]*tmpFx[87] + tmpQ2[124]*tmpFx[101] + tmpQ2[125]*tmpFx[115] + tmpQ2[126]*tmpFx[129] + tmpQ2[127]*tmpFx[143] + tmpQ2[128]*tmpFx[157] + tmpQ2[129]*tmpFx[171];
tmpQ1[130] = + tmpQ2[117]*tmpFx[4] + tmpQ2[118]*tmpFx[18] + tmpQ2[119]*tmpFx[32] + tmpQ2[120]*tmpFx[46] + tmpQ2[121]*tmpFx[60] + tmpQ2[122]*tmpFx[74] + tmpQ2[123]*tmpFx[88] + tmpQ2[124]*tmpFx[102] + tmpQ2[125]*tmpFx[116] + tmpQ2[126]*tmpFx[130] + tmpQ2[127]*tmpFx[144] + tmpQ2[128]*tmpFx[158] + tmpQ2[129]*tmpFx[172];
tmpQ1[131] = + tmpQ2[117]*tmpFx[5] + tmpQ2[118]*tmpFx[19] + tmpQ2[119]*tmpFx[33] + tmpQ2[120]*tmpFx[47] + tmpQ2[121]*tmpFx[61] + tmpQ2[122]*tmpFx[75] + tmpQ2[123]*tmpFx[89] + tmpQ2[124]*tmpFx[103] + tmpQ2[125]*tmpFx[117] + tmpQ2[126]*tmpFx[131] + tmpQ2[127]*tmpFx[145] + tmpQ2[128]*tmpFx[159] + tmpQ2[129]*tmpFx[173];
tmpQ1[132] = + tmpQ2[117]*tmpFx[6] + tmpQ2[118]*tmpFx[20] + tmpQ2[119]*tmpFx[34] + tmpQ2[120]*tmpFx[48] + tmpQ2[121]*tmpFx[62] + tmpQ2[122]*tmpFx[76] + tmpQ2[123]*tmpFx[90] + tmpQ2[124]*tmpFx[104] + tmpQ2[125]*tmpFx[118] + tmpQ2[126]*tmpFx[132] + tmpQ2[127]*tmpFx[146] + tmpQ2[128]*tmpFx[160] + tmpQ2[129]*tmpFx[174];
tmpQ1[133] = + tmpQ2[117]*tmpFx[7] + tmpQ2[118]*tmpFx[21] + tmpQ2[119]*tmpFx[35] + tmpQ2[120]*tmpFx[49] + tmpQ2[121]*tmpFx[63] + tmpQ2[122]*tmpFx[77] + tmpQ2[123]*tmpFx[91] + tmpQ2[124]*tmpFx[105] + tmpQ2[125]*tmpFx[119] + tmpQ2[126]*tmpFx[133] + tmpQ2[127]*tmpFx[147] + tmpQ2[128]*tmpFx[161] + tmpQ2[129]*tmpFx[175];
tmpQ1[134] = + tmpQ2[117]*tmpFx[8] + tmpQ2[118]*tmpFx[22] + tmpQ2[119]*tmpFx[36] + tmpQ2[120]*tmpFx[50] + tmpQ2[121]*tmpFx[64] + tmpQ2[122]*tmpFx[78] + tmpQ2[123]*tmpFx[92] + tmpQ2[124]*tmpFx[106] + tmpQ2[125]*tmpFx[120] + tmpQ2[126]*tmpFx[134] + tmpQ2[127]*tmpFx[148] + tmpQ2[128]*tmpFx[162] + tmpQ2[129]*tmpFx[176];
tmpQ1[135] = + tmpQ2[117]*tmpFx[9] + tmpQ2[118]*tmpFx[23] + tmpQ2[119]*tmpFx[37] + tmpQ2[120]*tmpFx[51] + tmpQ2[121]*tmpFx[65] + tmpQ2[122]*tmpFx[79] + tmpQ2[123]*tmpFx[93] + tmpQ2[124]*tmpFx[107] + tmpQ2[125]*tmpFx[121] + tmpQ2[126]*tmpFx[135] + tmpQ2[127]*tmpFx[149] + tmpQ2[128]*tmpFx[163] + tmpQ2[129]*tmpFx[177];
tmpQ1[136] = + tmpQ2[117]*tmpFx[10] + tmpQ2[118]*tmpFx[24] + tmpQ2[119]*tmpFx[38] + tmpQ2[120]*tmpFx[52] + tmpQ2[121]*tmpFx[66] + tmpQ2[122]*tmpFx[80] + tmpQ2[123]*tmpFx[94] + tmpQ2[124]*tmpFx[108] + tmpQ2[125]*tmpFx[122] + tmpQ2[126]*tmpFx[136] + tmpQ2[127]*tmpFx[150] + tmpQ2[128]*tmpFx[164] + tmpQ2[129]*tmpFx[178];
tmpQ1[137] = + tmpQ2[117]*tmpFx[11] + tmpQ2[118]*tmpFx[25] + tmpQ2[119]*tmpFx[39] + tmpQ2[120]*tmpFx[53] + tmpQ2[121]*tmpFx[67] + tmpQ2[122]*tmpFx[81] + tmpQ2[123]*tmpFx[95] + tmpQ2[124]*tmpFx[109] + tmpQ2[125]*tmpFx[123] + tmpQ2[126]*tmpFx[137] + tmpQ2[127]*tmpFx[151] + tmpQ2[128]*tmpFx[165] + tmpQ2[129]*tmpFx[179];
tmpQ1[138] = + tmpQ2[117]*tmpFx[12] + tmpQ2[118]*tmpFx[26] + tmpQ2[119]*tmpFx[40] + tmpQ2[120]*tmpFx[54] + tmpQ2[121]*tmpFx[68] + tmpQ2[122]*tmpFx[82] + tmpQ2[123]*tmpFx[96] + tmpQ2[124]*tmpFx[110] + tmpQ2[125]*tmpFx[124] + tmpQ2[126]*tmpFx[138] + tmpQ2[127]*tmpFx[152] + tmpQ2[128]*tmpFx[166] + tmpQ2[129]*tmpFx[180];
tmpQ1[139] = + tmpQ2[117]*tmpFx[13] + tmpQ2[118]*tmpFx[27] + tmpQ2[119]*tmpFx[41] + tmpQ2[120]*tmpFx[55] + tmpQ2[121]*tmpFx[69] + tmpQ2[122]*tmpFx[83] + tmpQ2[123]*tmpFx[97] + tmpQ2[124]*tmpFx[111] + tmpQ2[125]*tmpFx[125] + tmpQ2[126]*tmpFx[139] + tmpQ2[127]*tmpFx[153] + tmpQ2[128]*tmpFx[167] + tmpQ2[129]*tmpFx[181];
tmpQ1[140] = + tmpQ2[130]*tmpFx[0] + tmpQ2[131]*tmpFx[14] + tmpQ2[132]*tmpFx[28] + tmpQ2[133]*tmpFx[42] + tmpQ2[134]*tmpFx[56] + tmpQ2[135]*tmpFx[70] + tmpQ2[136]*tmpFx[84] + tmpQ2[137]*tmpFx[98] + tmpQ2[138]*tmpFx[112] + tmpQ2[139]*tmpFx[126] + tmpQ2[140]*tmpFx[140] + tmpQ2[141]*tmpFx[154] + tmpQ2[142]*tmpFx[168];
tmpQ1[141] = + tmpQ2[130]*tmpFx[1] + tmpQ2[131]*tmpFx[15] + tmpQ2[132]*tmpFx[29] + tmpQ2[133]*tmpFx[43] + tmpQ2[134]*tmpFx[57] + tmpQ2[135]*tmpFx[71] + tmpQ2[136]*tmpFx[85] + tmpQ2[137]*tmpFx[99] + tmpQ2[138]*tmpFx[113] + tmpQ2[139]*tmpFx[127] + tmpQ2[140]*tmpFx[141] + tmpQ2[141]*tmpFx[155] + tmpQ2[142]*tmpFx[169];
tmpQ1[142] = + tmpQ2[130]*tmpFx[2] + tmpQ2[131]*tmpFx[16] + tmpQ2[132]*tmpFx[30] + tmpQ2[133]*tmpFx[44] + tmpQ2[134]*tmpFx[58] + tmpQ2[135]*tmpFx[72] + tmpQ2[136]*tmpFx[86] + tmpQ2[137]*tmpFx[100] + tmpQ2[138]*tmpFx[114] + tmpQ2[139]*tmpFx[128] + tmpQ2[140]*tmpFx[142] + tmpQ2[141]*tmpFx[156] + tmpQ2[142]*tmpFx[170];
tmpQ1[143] = + tmpQ2[130]*tmpFx[3] + tmpQ2[131]*tmpFx[17] + tmpQ2[132]*tmpFx[31] + tmpQ2[133]*tmpFx[45] + tmpQ2[134]*tmpFx[59] + tmpQ2[135]*tmpFx[73] + tmpQ2[136]*tmpFx[87] + tmpQ2[137]*tmpFx[101] + tmpQ2[138]*tmpFx[115] + tmpQ2[139]*tmpFx[129] + tmpQ2[140]*tmpFx[143] + tmpQ2[141]*tmpFx[157] + tmpQ2[142]*tmpFx[171];
tmpQ1[144] = + tmpQ2[130]*tmpFx[4] + tmpQ2[131]*tmpFx[18] + tmpQ2[132]*tmpFx[32] + tmpQ2[133]*tmpFx[46] + tmpQ2[134]*tmpFx[60] + tmpQ2[135]*tmpFx[74] + tmpQ2[136]*tmpFx[88] + tmpQ2[137]*tmpFx[102] + tmpQ2[138]*tmpFx[116] + tmpQ2[139]*tmpFx[130] + tmpQ2[140]*tmpFx[144] + tmpQ2[141]*tmpFx[158] + tmpQ2[142]*tmpFx[172];
tmpQ1[145] = + tmpQ2[130]*tmpFx[5] + tmpQ2[131]*tmpFx[19] + tmpQ2[132]*tmpFx[33] + tmpQ2[133]*tmpFx[47] + tmpQ2[134]*tmpFx[61] + tmpQ2[135]*tmpFx[75] + tmpQ2[136]*tmpFx[89] + tmpQ2[137]*tmpFx[103] + tmpQ2[138]*tmpFx[117] + tmpQ2[139]*tmpFx[131] + tmpQ2[140]*tmpFx[145] + tmpQ2[141]*tmpFx[159] + tmpQ2[142]*tmpFx[173];
tmpQ1[146] = + tmpQ2[130]*tmpFx[6] + tmpQ2[131]*tmpFx[20] + tmpQ2[132]*tmpFx[34] + tmpQ2[133]*tmpFx[48] + tmpQ2[134]*tmpFx[62] + tmpQ2[135]*tmpFx[76] + tmpQ2[136]*tmpFx[90] + tmpQ2[137]*tmpFx[104] + tmpQ2[138]*tmpFx[118] + tmpQ2[139]*tmpFx[132] + tmpQ2[140]*tmpFx[146] + tmpQ2[141]*tmpFx[160] + tmpQ2[142]*tmpFx[174];
tmpQ1[147] = + tmpQ2[130]*tmpFx[7] + tmpQ2[131]*tmpFx[21] + tmpQ2[132]*tmpFx[35] + tmpQ2[133]*tmpFx[49] + tmpQ2[134]*tmpFx[63] + tmpQ2[135]*tmpFx[77] + tmpQ2[136]*tmpFx[91] + tmpQ2[137]*tmpFx[105] + tmpQ2[138]*tmpFx[119] + tmpQ2[139]*tmpFx[133] + tmpQ2[140]*tmpFx[147] + tmpQ2[141]*tmpFx[161] + tmpQ2[142]*tmpFx[175];
tmpQ1[148] = + tmpQ2[130]*tmpFx[8] + tmpQ2[131]*tmpFx[22] + tmpQ2[132]*tmpFx[36] + tmpQ2[133]*tmpFx[50] + tmpQ2[134]*tmpFx[64] + tmpQ2[135]*tmpFx[78] + tmpQ2[136]*tmpFx[92] + tmpQ2[137]*tmpFx[106] + tmpQ2[138]*tmpFx[120] + tmpQ2[139]*tmpFx[134] + tmpQ2[140]*tmpFx[148] + tmpQ2[141]*tmpFx[162] + tmpQ2[142]*tmpFx[176];
tmpQ1[149] = + tmpQ2[130]*tmpFx[9] + tmpQ2[131]*tmpFx[23] + tmpQ2[132]*tmpFx[37] + tmpQ2[133]*tmpFx[51] + tmpQ2[134]*tmpFx[65] + tmpQ2[135]*tmpFx[79] + tmpQ2[136]*tmpFx[93] + tmpQ2[137]*tmpFx[107] + tmpQ2[138]*tmpFx[121] + tmpQ2[139]*tmpFx[135] + tmpQ2[140]*tmpFx[149] + tmpQ2[141]*tmpFx[163] + tmpQ2[142]*tmpFx[177];
tmpQ1[150] = + tmpQ2[130]*tmpFx[10] + tmpQ2[131]*tmpFx[24] + tmpQ2[132]*tmpFx[38] + tmpQ2[133]*tmpFx[52] + tmpQ2[134]*tmpFx[66] + tmpQ2[135]*tmpFx[80] + tmpQ2[136]*tmpFx[94] + tmpQ2[137]*tmpFx[108] + tmpQ2[138]*tmpFx[122] + tmpQ2[139]*tmpFx[136] + tmpQ2[140]*tmpFx[150] + tmpQ2[141]*tmpFx[164] + tmpQ2[142]*tmpFx[178];
tmpQ1[151] = + tmpQ2[130]*tmpFx[11] + tmpQ2[131]*tmpFx[25] + tmpQ2[132]*tmpFx[39] + tmpQ2[133]*tmpFx[53] + tmpQ2[134]*tmpFx[67] + tmpQ2[135]*tmpFx[81] + tmpQ2[136]*tmpFx[95] + tmpQ2[137]*tmpFx[109] + tmpQ2[138]*tmpFx[123] + tmpQ2[139]*tmpFx[137] + tmpQ2[140]*tmpFx[151] + tmpQ2[141]*tmpFx[165] + tmpQ2[142]*tmpFx[179];
tmpQ1[152] = + tmpQ2[130]*tmpFx[12] + tmpQ2[131]*tmpFx[26] + tmpQ2[132]*tmpFx[40] + tmpQ2[133]*tmpFx[54] + tmpQ2[134]*tmpFx[68] + tmpQ2[135]*tmpFx[82] + tmpQ2[136]*tmpFx[96] + tmpQ2[137]*tmpFx[110] + tmpQ2[138]*tmpFx[124] + tmpQ2[139]*tmpFx[138] + tmpQ2[140]*tmpFx[152] + tmpQ2[141]*tmpFx[166] + tmpQ2[142]*tmpFx[180];
tmpQ1[153] = + tmpQ2[130]*tmpFx[13] + tmpQ2[131]*tmpFx[27] + tmpQ2[132]*tmpFx[41] + tmpQ2[133]*tmpFx[55] + tmpQ2[134]*tmpFx[69] + tmpQ2[135]*tmpFx[83] + tmpQ2[136]*tmpFx[97] + tmpQ2[137]*tmpFx[111] + tmpQ2[138]*tmpFx[125] + tmpQ2[139]*tmpFx[139] + tmpQ2[140]*tmpFx[153] + tmpQ2[141]*tmpFx[167] + tmpQ2[142]*tmpFx[181];
tmpQ1[154] = + tmpQ2[143]*tmpFx[0] + tmpQ2[144]*tmpFx[14] + tmpQ2[145]*tmpFx[28] + tmpQ2[146]*tmpFx[42] + tmpQ2[147]*tmpFx[56] + tmpQ2[148]*tmpFx[70] + tmpQ2[149]*tmpFx[84] + tmpQ2[150]*tmpFx[98] + tmpQ2[151]*tmpFx[112] + tmpQ2[152]*tmpFx[126] + tmpQ2[153]*tmpFx[140] + tmpQ2[154]*tmpFx[154] + tmpQ2[155]*tmpFx[168];
tmpQ1[155] = + tmpQ2[143]*tmpFx[1] + tmpQ2[144]*tmpFx[15] + tmpQ2[145]*tmpFx[29] + tmpQ2[146]*tmpFx[43] + tmpQ2[147]*tmpFx[57] + tmpQ2[148]*tmpFx[71] + tmpQ2[149]*tmpFx[85] + tmpQ2[150]*tmpFx[99] + tmpQ2[151]*tmpFx[113] + tmpQ2[152]*tmpFx[127] + tmpQ2[153]*tmpFx[141] + tmpQ2[154]*tmpFx[155] + tmpQ2[155]*tmpFx[169];
tmpQ1[156] = + tmpQ2[143]*tmpFx[2] + tmpQ2[144]*tmpFx[16] + tmpQ2[145]*tmpFx[30] + tmpQ2[146]*tmpFx[44] + tmpQ2[147]*tmpFx[58] + tmpQ2[148]*tmpFx[72] + tmpQ2[149]*tmpFx[86] + tmpQ2[150]*tmpFx[100] + tmpQ2[151]*tmpFx[114] + tmpQ2[152]*tmpFx[128] + tmpQ2[153]*tmpFx[142] + tmpQ2[154]*tmpFx[156] + tmpQ2[155]*tmpFx[170];
tmpQ1[157] = + tmpQ2[143]*tmpFx[3] + tmpQ2[144]*tmpFx[17] + tmpQ2[145]*tmpFx[31] + tmpQ2[146]*tmpFx[45] + tmpQ2[147]*tmpFx[59] + tmpQ2[148]*tmpFx[73] + tmpQ2[149]*tmpFx[87] + tmpQ2[150]*tmpFx[101] + tmpQ2[151]*tmpFx[115] + tmpQ2[152]*tmpFx[129] + tmpQ2[153]*tmpFx[143] + tmpQ2[154]*tmpFx[157] + tmpQ2[155]*tmpFx[171];
tmpQ1[158] = + tmpQ2[143]*tmpFx[4] + tmpQ2[144]*tmpFx[18] + tmpQ2[145]*tmpFx[32] + tmpQ2[146]*tmpFx[46] + tmpQ2[147]*tmpFx[60] + tmpQ2[148]*tmpFx[74] + tmpQ2[149]*tmpFx[88] + tmpQ2[150]*tmpFx[102] + tmpQ2[151]*tmpFx[116] + tmpQ2[152]*tmpFx[130] + tmpQ2[153]*tmpFx[144] + tmpQ2[154]*tmpFx[158] + tmpQ2[155]*tmpFx[172];
tmpQ1[159] = + tmpQ2[143]*tmpFx[5] + tmpQ2[144]*tmpFx[19] + tmpQ2[145]*tmpFx[33] + tmpQ2[146]*tmpFx[47] + tmpQ2[147]*tmpFx[61] + tmpQ2[148]*tmpFx[75] + tmpQ2[149]*tmpFx[89] + tmpQ2[150]*tmpFx[103] + tmpQ2[151]*tmpFx[117] + tmpQ2[152]*tmpFx[131] + tmpQ2[153]*tmpFx[145] + tmpQ2[154]*tmpFx[159] + tmpQ2[155]*tmpFx[173];
tmpQ1[160] = + tmpQ2[143]*tmpFx[6] + tmpQ2[144]*tmpFx[20] + tmpQ2[145]*tmpFx[34] + tmpQ2[146]*tmpFx[48] + tmpQ2[147]*tmpFx[62] + tmpQ2[148]*tmpFx[76] + tmpQ2[149]*tmpFx[90] + tmpQ2[150]*tmpFx[104] + tmpQ2[151]*tmpFx[118] + tmpQ2[152]*tmpFx[132] + tmpQ2[153]*tmpFx[146] + tmpQ2[154]*tmpFx[160] + tmpQ2[155]*tmpFx[174];
tmpQ1[161] = + tmpQ2[143]*tmpFx[7] + tmpQ2[144]*tmpFx[21] + tmpQ2[145]*tmpFx[35] + tmpQ2[146]*tmpFx[49] + tmpQ2[147]*tmpFx[63] + tmpQ2[148]*tmpFx[77] + tmpQ2[149]*tmpFx[91] + tmpQ2[150]*tmpFx[105] + tmpQ2[151]*tmpFx[119] + tmpQ2[152]*tmpFx[133] + tmpQ2[153]*tmpFx[147] + tmpQ2[154]*tmpFx[161] + tmpQ2[155]*tmpFx[175];
tmpQ1[162] = + tmpQ2[143]*tmpFx[8] + tmpQ2[144]*tmpFx[22] + tmpQ2[145]*tmpFx[36] + tmpQ2[146]*tmpFx[50] + tmpQ2[147]*tmpFx[64] + tmpQ2[148]*tmpFx[78] + tmpQ2[149]*tmpFx[92] + tmpQ2[150]*tmpFx[106] + tmpQ2[151]*tmpFx[120] + tmpQ2[152]*tmpFx[134] + tmpQ2[153]*tmpFx[148] + tmpQ2[154]*tmpFx[162] + tmpQ2[155]*tmpFx[176];
tmpQ1[163] = + tmpQ2[143]*tmpFx[9] + tmpQ2[144]*tmpFx[23] + tmpQ2[145]*tmpFx[37] + tmpQ2[146]*tmpFx[51] + tmpQ2[147]*tmpFx[65] + tmpQ2[148]*tmpFx[79] + tmpQ2[149]*tmpFx[93] + tmpQ2[150]*tmpFx[107] + tmpQ2[151]*tmpFx[121] + tmpQ2[152]*tmpFx[135] + tmpQ2[153]*tmpFx[149] + tmpQ2[154]*tmpFx[163] + tmpQ2[155]*tmpFx[177];
tmpQ1[164] = + tmpQ2[143]*tmpFx[10] + tmpQ2[144]*tmpFx[24] + tmpQ2[145]*tmpFx[38] + tmpQ2[146]*tmpFx[52] + tmpQ2[147]*tmpFx[66] + tmpQ2[148]*tmpFx[80] + tmpQ2[149]*tmpFx[94] + tmpQ2[150]*tmpFx[108] + tmpQ2[151]*tmpFx[122] + tmpQ2[152]*tmpFx[136] + tmpQ2[153]*tmpFx[150] + tmpQ2[154]*tmpFx[164] + tmpQ2[155]*tmpFx[178];
tmpQ1[165] = + tmpQ2[143]*tmpFx[11] + tmpQ2[144]*tmpFx[25] + tmpQ2[145]*tmpFx[39] + tmpQ2[146]*tmpFx[53] + tmpQ2[147]*tmpFx[67] + tmpQ2[148]*tmpFx[81] + tmpQ2[149]*tmpFx[95] + tmpQ2[150]*tmpFx[109] + tmpQ2[151]*tmpFx[123] + tmpQ2[152]*tmpFx[137] + tmpQ2[153]*tmpFx[151] + tmpQ2[154]*tmpFx[165] + tmpQ2[155]*tmpFx[179];
tmpQ1[166] = + tmpQ2[143]*tmpFx[12] + tmpQ2[144]*tmpFx[26] + tmpQ2[145]*tmpFx[40] + tmpQ2[146]*tmpFx[54] + tmpQ2[147]*tmpFx[68] + tmpQ2[148]*tmpFx[82] + tmpQ2[149]*tmpFx[96] + tmpQ2[150]*tmpFx[110] + tmpQ2[151]*tmpFx[124] + tmpQ2[152]*tmpFx[138] + tmpQ2[153]*tmpFx[152] + tmpQ2[154]*tmpFx[166] + tmpQ2[155]*tmpFx[180];
tmpQ1[167] = + tmpQ2[143]*tmpFx[13] + tmpQ2[144]*tmpFx[27] + tmpQ2[145]*tmpFx[41] + tmpQ2[146]*tmpFx[55] + tmpQ2[147]*tmpFx[69] + tmpQ2[148]*tmpFx[83] + tmpQ2[149]*tmpFx[97] + tmpQ2[150]*tmpFx[111] + tmpQ2[151]*tmpFx[125] + tmpQ2[152]*tmpFx[139] + tmpQ2[153]*tmpFx[153] + tmpQ2[154]*tmpFx[167] + tmpQ2[155]*tmpFx[181];
tmpQ1[168] = + tmpQ2[156]*tmpFx[0] + tmpQ2[157]*tmpFx[14] + tmpQ2[158]*tmpFx[28] + tmpQ2[159]*tmpFx[42] + tmpQ2[160]*tmpFx[56] + tmpQ2[161]*tmpFx[70] + tmpQ2[162]*tmpFx[84] + tmpQ2[163]*tmpFx[98] + tmpQ2[164]*tmpFx[112] + tmpQ2[165]*tmpFx[126] + tmpQ2[166]*tmpFx[140] + tmpQ2[167]*tmpFx[154] + tmpQ2[168]*tmpFx[168];
tmpQ1[169] = + tmpQ2[156]*tmpFx[1] + tmpQ2[157]*tmpFx[15] + tmpQ2[158]*tmpFx[29] + tmpQ2[159]*tmpFx[43] + tmpQ2[160]*tmpFx[57] + tmpQ2[161]*tmpFx[71] + tmpQ2[162]*tmpFx[85] + tmpQ2[163]*tmpFx[99] + tmpQ2[164]*tmpFx[113] + tmpQ2[165]*tmpFx[127] + tmpQ2[166]*tmpFx[141] + tmpQ2[167]*tmpFx[155] + tmpQ2[168]*tmpFx[169];
tmpQ1[170] = + tmpQ2[156]*tmpFx[2] + tmpQ2[157]*tmpFx[16] + tmpQ2[158]*tmpFx[30] + tmpQ2[159]*tmpFx[44] + tmpQ2[160]*tmpFx[58] + tmpQ2[161]*tmpFx[72] + tmpQ2[162]*tmpFx[86] + tmpQ2[163]*tmpFx[100] + tmpQ2[164]*tmpFx[114] + tmpQ2[165]*tmpFx[128] + tmpQ2[166]*tmpFx[142] + tmpQ2[167]*tmpFx[156] + tmpQ2[168]*tmpFx[170];
tmpQ1[171] = + tmpQ2[156]*tmpFx[3] + tmpQ2[157]*tmpFx[17] + tmpQ2[158]*tmpFx[31] + tmpQ2[159]*tmpFx[45] + tmpQ2[160]*tmpFx[59] + tmpQ2[161]*tmpFx[73] + tmpQ2[162]*tmpFx[87] + tmpQ2[163]*tmpFx[101] + tmpQ2[164]*tmpFx[115] + tmpQ2[165]*tmpFx[129] + tmpQ2[166]*tmpFx[143] + tmpQ2[167]*tmpFx[157] + tmpQ2[168]*tmpFx[171];
tmpQ1[172] = + tmpQ2[156]*tmpFx[4] + tmpQ2[157]*tmpFx[18] + tmpQ2[158]*tmpFx[32] + tmpQ2[159]*tmpFx[46] + tmpQ2[160]*tmpFx[60] + tmpQ2[161]*tmpFx[74] + tmpQ2[162]*tmpFx[88] + tmpQ2[163]*tmpFx[102] + tmpQ2[164]*tmpFx[116] + tmpQ2[165]*tmpFx[130] + tmpQ2[166]*tmpFx[144] + tmpQ2[167]*tmpFx[158] + tmpQ2[168]*tmpFx[172];
tmpQ1[173] = + tmpQ2[156]*tmpFx[5] + tmpQ2[157]*tmpFx[19] + tmpQ2[158]*tmpFx[33] + tmpQ2[159]*tmpFx[47] + tmpQ2[160]*tmpFx[61] + tmpQ2[161]*tmpFx[75] + tmpQ2[162]*tmpFx[89] + tmpQ2[163]*tmpFx[103] + tmpQ2[164]*tmpFx[117] + tmpQ2[165]*tmpFx[131] + tmpQ2[166]*tmpFx[145] + tmpQ2[167]*tmpFx[159] + tmpQ2[168]*tmpFx[173];
tmpQ1[174] = + tmpQ2[156]*tmpFx[6] + tmpQ2[157]*tmpFx[20] + tmpQ2[158]*tmpFx[34] + tmpQ2[159]*tmpFx[48] + tmpQ2[160]*tmpFx[62] + tmpQ2[161]*tmpFx[76] + tmpQ2[162]*tmpFx[90] + tmpQ2[163]*tmpFx[104] + tmpQ2[164]*tmpFx[118] + tmpQ2[165]*tmpFx[132] + tmpQ2[166]*tmpFx[146] + tmpQ2[167]*tmpFx[160] + tmpQ2[168]*tmpFx[174];
tmpQ1[175] = + tmpQ2[156]*tmpFx[7] + tmpQ2[157]*tmpFx[21] + tmpQ2[158]*tmpFx[35] + tmpQ2[159]*tmpFx[49] + tmpQ2[160]*tmpFx[63] + tmpQ2[161]*tmpFx[77] + tmpQ2[162]*tmpFx[91] + tmpQ2[163]*tmpFx[105] + tmpQ2[164]*tmpFx[119] + tmpQ2[165]*tmpFx[133] + tmpQ2[166]*tmpFx[147] + tmpQ2[167]*tmpFx[161] + tmpQ2[168]*tmpFx[175];
tmpQ1[176] = + tmpQ2[156]*tmpFx[8] + tmpQ2[157]*tmpFx[22] + tmpQ2[158]*tmpFx[36] + tmpQ2[159]*tmpFx[50] + tmpQ2[160]*tmpFx[64] + tmpQ2[161]*tmpFx[78] + tmpQ2[162]*tmpFx[92] + tmpQ2[163]*tmpFx[106] + tmpQ2[164]*tmpFx[120] + tmpQ2[165]*tmpFx[134] + tmpQ2[166]*tmpFx[148] + tmpQ2[167]*tmpFx[162] + tmpQ2[168]*tmpFx[176];
tmpQ1[177] = + tmpQ2[156]*tmpFx[9] + tmpQ2[157]*tmpFx[23] + tmpQ2[158]*tmpFx[37] + tmpQ2[159]*tmpFx[51] + tmpQ2[160]*tmpFx[65] + tmpQ2[161]*tmpFx[79] + tmpQ2[162]*tmpFx[93] + tmpQ2[163]*tmpFx[107] + tmpQ2[164]*tmpFx[121] + tmpQ2[165]*tmpFx[135] + tmpQ2[166]*tmpFx[149] + tmpQ2[167]*tmpFx[163] + tmpQ2[168]*tmpFx[177];
tmpQ1[178] = + tmpQ2[156]*tmpFx[10] + tmpQ2[157]*tmpFx[24] + tmpQ2[158]*tmpFx[38] + tmpQ2[159]*tmpFx[52] + tmpQ2[160]*tmpFx[66] + tmpQ2[161]*tmpFx[80] + tmpQ2[162]*tmpFx[94] + tmpQ2[163]*tmpFx[108] + tmpQ2[164]*tmpFx[122] + tmpQ2[165]*tmpFx[136] + tmpQ2[166]*tmpFx[150] + tmpQ2[167]*tmpFx[164] + tmpQ2[168]*tmpFx[178];
tmpQ1[179] = + tmpQ2[156]*tmpFx[11] + tmpQ2[157]*tmpFx[25] + tmpQ2[158]*tmpFx[39] + tmpQ2[159]*tmpFx[53] + tmpQ2[160]*tmpFx[67] + tmpQ2[161]*tmpFx[81] + tmpQ2[162]*tmpFx[95] + tmpQ2[163]*tmpFx[109] + tmpQ2[164]*tmpFx[123] + tmpQ2[165]*tmpFx[137] + tmpQ2[166]*tmpFx[151] + tmpQ2[167]*tmpFx[165] + tmpQ2[168]*tmpFx[179];
tmpQ1[180] = + tmpQ2[156]*tmpFx[12] + tmpQ2[157]*tmpFx[26] + tmpQ2[158]*tmpFx[40] + tmpQ2[159]*tmpFx[54] + tmpQ2[160]*tmpFx[68] + tmpQ2[161]*tmpFx[82] + tmpQ2[162]*tmpFx[96] + tmpQ2[163]*tmpFx[110] + tmpQ2[164]*tmpFx[124] + tmpQ2[165]*tmpFx[138] + tmpQ2[166]*tmpFx[152] + tmpQ2[167]*tmpFx[166] + tmpQ2[168]*tmpFx[180];
tmpQ1[181] = + tmpQ2[156]*tmpFx[13] + tmpQ2[157]*tmpFx[27] + tmpQ2[158]*tmpFx[41] + tmpQ2[159]*tmpFx[55] + tmpQ2[160]*tmpFx[69] + tmpQ2[161]*tmpFx[83] + tmpQ2[162]*tmpFx[97] + tmpQ2[163]*tmpFx[111] + tmpQ2[164]*tmpFx[125] + tmpQ2[165]*tmpFx[139] + tmpQ2[166]*tmpFx[153] + tmpQ2[167]*tmpFx[167] + tmpQ2[168]*tmpFx[181];
tmpQ1[182] = + tmpQ2[169]*tmpFx[0] + tmpQ2[170]*tmpFx[14] + tmpQ2[171]*tmpFx[28] + tmpQ2[172]*tmpFx[42] + tmpQ2[173]*tmpFx[56] + tmpQ2[174]*tmpFx[70] + tmpQ2[175]*tmpFx[84] + tmpQ2[176]*tmpFx[98] + tmpQ2[177]*tmpFx[112] + tmpQ2[178]*tmpFx[126] + tmpQ2[179]*tmpFx[140] + tmpQ2[180]*tmpFx[154] + tmpQ2[181]*tmpFx[168];
tmpQ1[183] = + tmpQ2[169]*tmpFx[1] + tmpQ2[170]*tmpFx[15] + tmpQ2[171]*tmpFx[29] + tmpQ2[172]*tmpFx[43] + tmpQ2[173]*tmpFx[57] + tmpQ2[174]*tmpFx[71] + tmpQ2[175]*tmpFx[85] + tmpQ2[176]*tmpFx[99] + tmpQ2[177]*tmpFx[113] + tmpQ2[178]*tmpFx[127] + tmpQ2[179]*tmpFx[141] + tmpQ2[180]*tmpFx[155] + tmpQ2[181]*tmpFx[169];
tmpQ1[184] = + tmpQ2[169]*tmpFx[2] + tmpQ2[170]*tmpFx[16] + tmpQ2[171]*tmpFx[30] + tmpQ2[172]*tmpFx[44] + tmpQ2[173]*tmpFx[58] + tmpQ2[174]*tmpFx[72] + tmpQ2[175]*tmpFx[86] + tmpQ2[176]*tmpFx[100] + tmpQ2[177]*tmpFx[114] + tmpQ2[178]*tmpFx[128] + tmpQ2[179]*tmpFx[142] + tmpQ2[180]*tmpFx[156] + tmpQ2[181]*tmpFx[170];
tmpQ1[185] = + tmpQ2[169]*tmpFx[3] + tmpQ2[170]*tmpFx[17] + tmpQ2[171]*tmpFx[31] + tmpQ2[172]*tmpFx[45] + tmpQ2[173]*tmpFx[59] + tmpQ2[174]*tmpFx[73] + tmpQ2[175]*tmpFx[87] + tmpQ2[176]*tmpFx[101] + tmpQ2[177]*tmpFx[115] + tmpQ2[178]*tmpFx[129] + tmpQ2[179]*tmpFx[143] + tmpQ2[180]*tmpFx[157] + tmpQ2[181]*tmpFx[171];
tmpQ1[186] = + tmpQ2[169]*tmpFx[4] + tmpQ2[170]*tmpFx[18] + tmpQ2[171]*tmpFx[32] + tmpQ2[172]*tmpFx[46] + tmpQ2[173]*tmpFx[60] + tmpQ2[174]*tmpFx[74] + tmpQ2[175]*tmpFx[88] + tmpQ2[176]*tmpFx[102] + tmpQ2[177]*tmpFx[116] + tmpQ2[178]*tmpFx[130] + tmpQ2[179]*tmpFx[144] + tmpQ2[180]*tmpFx[158] + tmpQ2[181]*tmpFx[172];
tmpQ1[187] = + tmpQ2[169]*tmpFx[5] + tmpQ2[170]*tmpFx[19] + tmpQ2[171]*tmpFx[33] + tmpQ2[172]*tmpFx[47] + tmpQ2[173]*tmpFx[61] + tmpQ2[174]*tmpFx[75] + tmpQ2[175]*tmpFx[89] + tmpQ2[176]*tmpFx[103] + tmpQ2[177]*tmpFx[117] + tmpQ2[178]*tmpFx[131] + tmpQ2[179]*tmpFx[145] + tmpQ2[180]*tmpFx[159] + tmpQ2[181]*tmpFx[173];
tmpQ1[188] = + tmpQ2[169]*tmpFx[6] + tmpQ2[170]*tmpFx[20] + tmpQ2[171]*tmpFx[34] + tmpQ2[172]*tmpFx[48] + tmpQ2[173]*tmpFx[62] + tmpQ2[174]*tmpFx[76] + tmpQ2[175]*tmpFx[90] + tmpQ2[176]*tmpFx[104] + tmpQ2[177]*tmpFx[118] + tmpQ2[178]*tmpFx[132] + tmpQ2[179]*tmpFx[146] + tmpQ2[180]*tmpFx[160] + tmpQ2[181]*tmpFx[174];
tmpQ1[189] = + tmpQ2[169]*tmpFx[7] + tmpQ2[170]*tmpFx[21] + tmpQ2[171]*tmpFx[35] + tmpQ2[172]*tmpFx[49] + tmpQ2[173]*tmpFx[63] + tmpQ2[174]*tmpFx[77] + tmpQ2[175]*tmpFx[91] + tmpQ2[176]*tmpFx[105] + tmpQ2[177]*tmpFx[119] + tmpQ2[178]*tmpFx[133] + tmpQ2[179]*tmpFx[147] + tmpQ2[180]*tmpFx[161] + tmpQ2[181]*tmpFx[175];
tmpQ1[190] = + tmpQ2[169]*tmpFx[8] + tmpQ2[170]*tmpFx[22] + tmpQ2[171]*tmpFx[36] + tmpQ2[172]*tmpFx[50] + tmpQ2[173]*tmpFx[64] + tmpQ2[174]*tmpFx[78] + tmpQ2[175]*tmpFx[92] + tmpQ2[176]*tmpFx[106] + tmpQ2[177]*tmpFx[120] + tmpQ2[178]*tmpFx[134] + tmpQ2[179]*tmpFx[148] + tmpQ2[180]*tmpFx[162] + tmpQ2[181]*tmpFx[176];
tmpQ1[191] = + tmpQ2[169]*tmpFx[9] + tmpQ2[170]*tmpFx[23] + tmpQ2[171]*tmpFx[37] + tmpQ2[172]*tmpFx[51] + tmpQ2[173]*tmpFx[65] + tmpQ2[174]*tmpFx[79] + tmpQ2[175]*tmpFx[93] + tmpQ2[176]*tmpFx[107] + tmpQ2[177]*tmpFx[121] + tmpQ2[178]*tmpFx[135] + tmpQ2[179]*tmpFx[149] + tmpQ2[180]*tmpFx[163] + tmpQ2[181]*tmpFx[177];
tmpQ1[192] = + tmpQ2[169]*tmpFx[10] + tmpQ2[170]*tmpFx[24] + tmpQ2[171]*tmpFx[38] + tmpQ2[172]*tmpFx[52] + tmpQ2[173]*tmpFx[66] + tmpQ2[174]*tmpFx[80] + tmpQ2[175]*tmpFx[94] + tmpQ2[176]*tmpFx[108] + tmpQ2[177]*tmpFx[122] + tmpQ2[178]*tmpFx[136] + tmpQ2[179]*tmpFx[150] + tmpQ2[180]*tmpFx[164] + tmpQ2[181]*tmpFx[178];
tmpQ1[193] = + tmpQ2[169]*tmpFx[11] + tmpQ2[170]*tmpFx[25] + tmpQ2[171]*tmpFx[39] + tmpQ2[172]*tmpFx[53] + tmpQ2[173]*tmpFx[67] + tmpQ2[174]*tmpFx[81] + tmpQ2[175]*tmpFx[95] + tmpQ2[176]*tmpFx[109] + tmpQ2[177]*tmpFx[123] + tmpQ2[178]*tmpFx[137] + tmpQ2[179]*tmpFx[151] + tmpQ2[180]*tmpFx[165] + tmpQ2[181]*tmpFx[179];
tmpQ1[194] = + tmpQ2[169]*tmpFx[12] + tmpQ2[170]*tmpFx[26] + tmpQ2[171]*tmpFx[40] + tmpQ2[172]*tmpFx[54] + tmpQ2[173]*tmpFx[68] + tmpQ2[174]*tmpFx[82] + tmpQ2[175]*tmpFx[96] + tmpQ2[176]*tmpFx[110] + tmpQ2[177]*tmpFx[124] + tmpQ2[178]*tmpFx[138] + tmpQ2[179]*tmpFx[152] + tmpQ2[180]*tmpFx[166] + tmpQ2[181]*tmpFx[180];
tmpQ1[195] = + tmpQ2[169]*tmpFx[13] + tmpQ2[170]*tmpFx[27] + tmpQ2[171]*tmpFx[41] + tmpQ2[172]*tmpFx[55] + tmpQ2[173]*tmpFx[69] + tmpQ2[174]*tmpFx[83] + tmpQ2[175]*tmpFx[97] + tmpQ2[176]*tmpFx[111] + tmpQ2[177]*tmpFx[125] + tmpQ2[178]*tmpFx[139] + tmpQ2[179]*tmpFx[153] + tmpQ2[180]*tmpFx[167] + tmpQ2[181]*tmpFx[181];
}

void nmpc_setObjR1R2( real_t* const tmpObjS, real_t* const tmpR1, real_t* const tmpR2 )
{
tmpR2[0] = +tmpObjS[117];
tmpR2[1] = +tmpObjS[118];
tmpR2[2] = +tmpObjS[119];
tmpR2[3] = +tmpObjS[120];
tmpR2[4] = +tmpObjS[121];
tmpR2[5] = +tmpObjS[122];
tmpR2[6] = +tmpObjS[123];
tmpR2[7] = +tmpObjS[124];
tmpR2[8] = +tmpObjS[125];
tmpR2[9] = +tmpObjS[126];
tmpR2[10] = +tmpObjS[127];
tmpR2[11] = +tmpObjS[128];
tmpR2[12] = +tmpObjS[129];
tmpR2[13] = +tmpObjS[130];
tmpR2[14] = +tmpObjS[131];
tmpR2[15] = +tmpObjS[132];
tmpR2[16] = +tmpObjS[133];
tmpR2[17] = +tmpObjS[134];
tmpR2[18] = +tmpObjS[135];
tmpR2[19] = +tmpObjS[136];
tmpR2[20] = +tmpObjS[137];
tmpR2[21] = +tmpObjS[138];
tmpR2[22] = +tmpObjS[139];
tmpR2[23] = +tmpObjS[140];
tmpR2[24] = +tmpObjS[141];
tmpR2[25] = +tmpObjS[142];
tmpR2[26] = +tmpObjS[143];
tmpR2[27] = +tmpObjS[144];
tmpR2[28] = +tmpObjS[145];
tmpR2[29] = +tmpObjS[146];
tmpR2[30] = +tmpObjS[147];
tmpR2[31] = +tmpObjS[148];
tmpR2[32] = +tmpObjS[149];
tmpR2[33] = +tmpObjS[150];
tmpR2[34] = +tmpObjS[151];
tmpR2[35] = +tmpObjS[152];
tmpR2[36] = +tmpObjS[153];
tmpR2[37] = +tmpObjS[154];
tmpR2[38] = +tmpObjS[155];
tmpR2[39] = +tmpObjS[156];
tmpR2[40] = +tmpObjS[157];
tmpR2[41] = +tmpObjS[158];
tmpR2[42] = +tmpObjS[159];
tmpR2[43] = +tmpObjS[160];
tmpR2[44] = +tmpObjS[161];
tmpR2[45] = +tmpObjS[162];
tmpR2[46] = +tmpObjS[163];
tmpR2[47] = +tmpObjS[164];
tmpR2[48] = +tmpObjS[165];
tmpR2[49] = +tmpObjS[166];
tmpR2[50] = +tmpObjS[167];
tmpR2[51] = +tmpObjS[168];
tmpR1[0] = + tmpR2[9];
tmpR1[1] = + tmpR2[10];
tmpR1[2] = + tmpR2[11];
tmpR1[3] = + tmpR2[12];
tmpR1[4] = + tmpR2[22];
tmpR1[5] = + tmpR2[23];
tmpR1[6] = + tmpR2[24];
tmpR1[7] = + tmpR2[25];
tmpR1[8] = + tmpR2[35];
tmpR1[9] = + tmpR2[36];
tmpR1[10] = + tmpR2[37];
tmpR1[11] = + tmpR2[38];
tmpR1[12] = + tmpR2[48];
tmpR1[13] = + tmpR2[49];
tmpR1[14] = + tmpR2[50];
tmpR1[15] = + tmpR2[51];
}

void nmpc_setObjQN1QN2( real_t* const tmpObjSEndTerm, real_t* const tmpQN1, real_t* const tmpQN2 )
{
tmpQN2[0] = +tmpObjSEndTerm[0];
tmpQN2[1] = +tmpObjSEndTerm[1];
tmpQN2[2] = +tmpObjSEndTerm[2];
tmpQN2[3] = +tmpObjSEndTerm[3];
tmpQN2[4] = +tmpObjSEndTerm[4];
tmpQN2[5] = +tmpObjSEndTerm[5];
tmpQN2[6] = +tmpObjSEndTerm[6];
tmpQN2[7] = +tmpObjSEndTerm[7];
tmpQN2[8] = +tmpObjSEndTerm[8];
tmpQN2[9] = +tmpObjSEndTerm[9];
tmpQN2[10] = +tmpObjSEndTerm[10];
tmpQN2[11] = +tmpObjSEndTerm[11];
tmpQN2[12] = +tmpObjSEndTerm[12];
tmpQN2[13] = +tmpObjSEndTerm[13];
tmpQN2[14] = +tmpObjSEndTerm[14];
tmpQN2[15] = +tmpObjSEndTerm[15];
tmpQN2[16] = +tmpObjSEndTerm[16];
tmpQN2[17] = +tmpObjSEndTerm[17];
tmpQN2[18] = +tmpObjSEndTerm[18];
tmpQN2[19] = +tmpObjSEndTerm[19];
tmpQN2[20] = +tmpObjSEndTerm[20];
tmpQN2[21] = +tmpObjSEndTerm[21];
tmpQN2[22] = +tmpObjSEndTerm[22];
tmpQN2[23] = +tmpObjSEndTerm[23];
tmpQN2[24] = +tmpObjSEndTerm[24];
tmpQN2[25] = +tmpObjSEndTerm[25];
tmpQN2[26] = +tmpObjSEndTerm[26];
tmpQN2[27] = +tmpObjSEndTerm[27];
tmpQN2[28] = +tmpObjSEndTerm[28];
tmpQN2[29] = +tmpObjSEndTerm[29];
tmpQN2[30] = +tmpObjSEndTerm[30];
tmpQN2[31] = +tmpObjSEndTerm[31];
tmpQN2[32] = +tmpObjSEndTerm[32];
tmpQN2[33] = +tmpObjSEndTerm[33];
tmpQN2[34] = +tmpObjSEndTerm[34];
tmpQN2[35] = +tmpObjSEndTerm[35];
tmpQN2[36] = +tmpObjSEndTerm[36];
tmpQN2[37] = +tmpObjSEndTerm[37];
tmpQN2[38] = +tmpObjSEndTerm[38];
tmpQN2[39] = +tmpObjSEndTerm[39];
tmpQN2[40] = +tmpObjSEndTerm[40];
tmpQN2[41] = +tmpObjSEndTerm[41];
tmpQN2[42] = +tmpObjSEndTerm[42];
tmpQN2[43] = +tmpObjSEndTerm[43];
tmpQN2[44] = +tmpObjSEndTerm[44];
tmpQN2[45] = +tmpObjSEndTerm[45];
tmpQN2[46] = +tmpObjSEndTerm[46];
tmpQN2[47] = +tmpObjSEndTerm[47];
tmpQN2[48] = +tmpObjSEndTerm[48];
tmpQN2[49] = +tmpObjSEndTerm[49];
tmpQN2[50] = +tmpObjSEndTerm[50];
tmpQN2[51] = +tmpObjSEndTerm[51];
tmpQN2[52] = +tmpObjSEndTerm[52];
tmpQN2[53] = +tmpObjSEndTerm[53];
tmpQN2[54] = +tmpObjSEndTerm[54];
tmpQN2[55] = +tmpObjSEndTerm[55];
tmpQN2[56] = +tmpObjSEndTerm[56];
tmpQN2[57] = +tmpObjSEndTerm[57];
tmpQN2[58] = +tmpObjSEndTerm[58];
tmpQN2[59] = +tmpObjSEndTerm[59];
tmpQN2[60] = +tmpObjSEndTerm[60];
tmpQN2[61] = +tmpObjSEndTerm[61];
tmpQN2[62] = +tmpObjSEndTerm[62];
tmpQN2[63] = +tmpObjSEndTerm[63];
tmpQN2[64] = 0.0;
;
tmpQN2[65] = 0.0;
;
tmpQN2[66] = 0.0;
;
tmpQN2[67] = 0.0;
;
tmpQN2[68] = 0.0;
;
tmpQN2[69] = 0.0;
;
tmpQN2[70] = 0.0;
;
tmpQN2[71] = 0.0;
;
tmpQN2[72] = 0.0;
;
tmpQN2[73] = 0.0;
;
tmpQN2[74] = 0.0;
;
tmpQN2[75] = 0.0;
;
tmpQN2[76] = 0.0;
;
tmpQN2[77] = 0.0;
;
tmpQN2[78] = 0.0;
;
tmpQN2[79] = 0.0;
;
tmpQN2[80] = 0.0;
;
tmpQN2[81] = 0.0;
;
tmpQN2[82] = 0.0;
;
tmpQN2[83] = 0.0;
;
tmpQN2[84] = 0.0;
;
tmpQN2[85] = 0.0;
;
tmpQN2[86] = 0.0;
;
tmpQN2[87] = 0.0;
;
tmpQN2[88] = 0.0;
;
tmpQN2[89] = 0.0;
;
tmpQN2[90] = 0.0;
;
tmpQN2[91] = 0.0;
;
tmpQN2[92] = 0.0;
;
tmpQN2[93] = 0.0;
;
tmpQN2[94] = 0.0;
;
tmpQN2[95] = 0.0;
;
tmpQN2[96] = 0.0;
;
tmpQN2[97] = 0.0;
;
tmpQN2[98] = 0.0;
;
tmpQN2[99] = 0.0;
;
tmpQN2[100] = 0.0;
;
tmpQN2[101] = 0.0;
;
tmpQN2[102] = 0.0;
;
tmpQN2[103] = 0.0;
;
tmpQN2[104] = 0.0;
;
tmpQN2[105] = 0.0;
;
tmpQN2[106] = 0.0;
;
tmpQN2[107] = 0.0;
;
tmpQN2[108] = 0.0;
;
tmpQN2[109] = 0.0;
;
tmpQN2[110] = 0.0;
;
tmpQN2[111] = 0.0;
;
tmpQN1[0] = + tmpQN2[0];
tmpQN1[1] = + tmpQN2[1];
tmpQN1[2] = + tmpQN2[2];
tmpQN1[3] = + tmpQN2[3];
tmpQN1[4] = + tmpQN2[4];
tmpQN1[5] = + tmpQN2[5];
tmpQN1[6] = + tmpQN2[6];
tmpQN1[7] = + tmpQN2[7];
tmpQN1[8] = 0.0;
;
tmpQN1[9] = 0.0;
;
tmpQN1[10] = 0.0;
;
tmpQN1[11] = 0.0;
;
tmpQN1[12] = 0.0;
;
tmpQN1[13] = 0.0;
;
tmpQN1[14] = + tmpQN2[8];
tmpQN1[15] = + tmpQN2[9];
tmpQN1[16] = + tmpQN2[10];
tmpQN1[17] = + tmpQN2[11];
tmpQN1[18] = + tmpQN2[12];
tmpQN1[19] = + tmpQN2[13];
tmpQN1[20] = + tmpQN2[14];
tmpQN1[21] = + tmpQN2[15];
tmpQN1[22] = 0.0;
;
tmpQN1[23] = 0.0;
;
tmpQN1[24] = 0.0;
;
tmpQN1[25] = 0.0;
;
tmpQN1[26] = 0.0;
;
tmpQN1[27] = 0.0;
;
tmpQN1[28] = + tmpQN2[16];
tmpQN1[29] = + tmpQN2[17];
tmpQN1[30] = + tmpQN2[18];
tmpQN1[31] = + tmpQN2[19];
tmpQN1[32] = + tmpQN2[20];
tmpQN1[33] = + tmpQN2[21];
tmpQN1[34] = + tmpQN2[22];
tmpQN1[35] = + tmpQN2[23];
tmpQN1[36] = 0.0;
;
tmpQN1[37] = 0.0;
;
tmpQN1[38] = 0.0;
;
tmpQN1[39] = 0.0;
;
tmpQN1[40] = 0.0;
;
tmpQN1[41] = 0.0;
;
tmpQN1[42] = + tmpQN2[24];
tmpQN1[43] = + tmpQN2[25];
tmpQN1[44] = + tmpQN2[26];
tmpQN1[45] = + tmpQN2[27];
tmpQN1[46] = + tmpQN2[28];
tmpQN1[47] = + tmpQN2[29];
tmpQN1[48] = + tmpQN2[30];
tmpQN1[49] = + tmpQN2[31];
tmpQN1[50] = 0.0;
;
tmpQN1[51] = 0.0;
;
tmpQN1[52] = 0.0;
;
tmpQN1[53] = 0.0;
;
tmpQN1[54] = 0.0;
;
tmpQN1[55] = 0.0;
;
tmpQN1[56] = + tmpQN2[32];
tmpQN1[57] = + tmpQN2[33];
tmpQN1[58] = + tmpQN2[34];
tmpQN1[59] = + tmpQN2[35];
tmpQN1[60] = + tmpQN2[36];
tmpQN1[61] = + tmpQN2[37];
tmpQN1[62] = + tmpQN2[38];
tmpQN1[63] = + tmpQN2[39];
tmpQN1[64] = 0.0;
;
tmpQN1[65] = 0.0;
;
tmpQN1[66] = 0.0;
;
tmpQN1[67] = 0.0;
;
tmpQN1[68] = 0.0;
;
tmpQN1[69] = 0.0;
;
tmpQN1[70] = + tmpQN2[40];
tmpQN1[71] = + tmpQN2[41];
tmpQN1[72] = + tmpQN2[42];
tmpQN1[73] = + tmpQN2[43];
tmpQN1[74] = + tmpQN2[44];
tmpQN1[75] = + tmpQN2[45];
tmpQN1[76] = + tmpQN2[46];
tmpQN1[77] = + tmpQN2[47];
tmpQN1[78] = 0.0;
;
tmpQN1[79] = 0.0;
;
tmpQN1[80] = 0.0;
;
tmpQN1[81] = 0.0;
;
tmpQN1[82] = 0.0;
;
tmpQN1[83] = 0.0;
;
tmpQN1[84] = + tmpQN2[48];
tmpQN1[85] = + tmpQN2[49];
tmpQN1[86] = + tmpQN2[50];
tmpQN1[87] = + tmpQN2[51];
tmpQN1[88] = + tmpQN2[52];
tmpQN1[89] = + tmpQN2[53];
tmpQN1[90] = + tmpQN2[54];
tmpQN1[91] = + tmpQN2[55];
tmpQN1[92] = 0.0;
;
tmpQN1[93] = 0.0;
;
tmpQN1[94] = 0.0;
;
tmpQN1[95] = 0.0;
;
tmpQN1[96] = 0.0;
;
tmpQN1[97] = 0.0;
;
tmpQN1[98] = + tmpQN2[56];
tmpQN1[99] = + tmpQN2[57];
tmpQN1[100] = + tmpQN2[58];
tmpQN1[101] = + tmpQN2[59];
tmpQN1[102] = + tmpQN2[60];
tmpQN1[103] = + tmpQN2[61];
tmpQN1[104] = + tmpQN2[62];
tmpQN1[105] = + tmpQN2[63];
tmpQN1[106] = 0.0;
;
tmpQN1[107] = 0.0;
;
tmpQN1[108] = 0.0;
;
tmpQN1[109] = 0.0;
;
tmpQN1[110] = 0.0;
;
tmpQN1[111] = 0.0;
;
tmpQN1[112] = + tmpQN2[64];
tmpQN1[113] = + tmpQN2[65];
tmpQN1[114] = + tmpQN2[66];
tmpQN1[115] = + tmpQN2[67];
tmpQN1[116] = + tmpQN2[68];
tmpQN1[117] = + tmpQN2[69];
tmpQN1[118] = + tmpQN2[70];
tmpQN1[119] = + tmpQN2[71];
tmpQN1[120] = 0.0;
;
tmpQN1[121] = 0.0;
;
tmpQN1[122] = 0.0;
;
tmpQN1[123] = 0.0;
;
tmpQN1[124] = 0.0;
;
tmpQN1[125] = 0.0;
;
tmpQN1[126] = + tmpQN2[72];
tmpQN1[127] = + tmpQN2[73];
tmpQN1[128] = + tmpQN2[74];
tmpQN1[129] = + tmpQN2[75];
tmpQN1[130] = + tmpQN2[76];
tmpQN1[131] = + tmpQN2[77];
tmpQN1[132] = + tmpQN2[78];
tmpQN1[133] = + tmpQN2[79];
tmpQN1[134] = 0.0;
;
tmpQN1[135] = 0.0;
;
tmpQN1[136] = 0.0;
;
tmpQN1[137] = 0.0;
;
tmpQN1[138] = 0.0;
;
tmpQN1[139] = 0.0;
;
tmpQN1[140] = + tmpQN2[80];
tmpQN1[141] = + tmpQN2[81];
tmpQN1[142] = + tmpQN2[82];
tmpQN1[143] = + tmpQN2[83];
tmpQN1[144] = + tmpQN2[84];
tmpQN1[145] = + tmpQN2[85];
tmpQN1[146] = + tmpQN2[86];
tmpQN1[147] = + tmpQN2[87];
tmpQN1[148] = 0.0;
;
tmpQN1[149] = 0.0;
;
tmpQN1[150] = 0.0;
;
tmpQN1[151] = 0.0;
;
tmpQN1[152] = 0.0;
;
tmpQN1[153] = 0.0;
;
tmpQN1[154] = + tmpQN2[88];
tmpQN1[155] = + tmpQN2[89];
tmpQN1[156] = + tmpQN2[90];
tmpQN1[157] = + tmpQN2[91];
tmpQN1[158] = + tmpQN2[92];
tmpQN1[159] = + tmpQN2[93];
tmpQN1[160] = + tmpQN2[94];
tmpQN1[161] = + tmpQN2[95];
tmpQN1[162] = 0.0;
;
tmpQN1[163] = 0.0;
;
tmpQN1[164] = 0.0;
;
tmpQN1[165] = 0.0;
;
tmpQN1[166] = 0.0;
;
tmpQN1[167] = 0.0;
;
tmpQN1[168] = + tmpQN2[96];
tmpQN1[169] = + tmpQN2[97];
tmpQN1[170] = + tmpQN2[98];
tmpQN1[171] = + tmpQN2[99];
tmpQN1[172] = + tmpQN2[100];
tmpQN1[173] = + tmpQN2[101];
tmpQN1[174] = + tmpQN2[102];
tmpQN1[175] = + tmpQN2[103];
tmpQN1[176] = 0.0;
;
tmpQN1[177] = 0.0;
;
tmpQN1[178] = 0.0;
;
tmpQN1[179] = 0.0;
;
tmpQN1[180] = 0.0;
;
tmpQN1[181] = 0.0;
;
tmpQN1[182] = + tmpQN2[104];
tmpQN1[183] = + tmpQN2[105];
tmpQN1[184] = + tmpQN2[106];
tmpQN1[185] = + tmpQN2[107];
tmpQN1[186] = + tmpQN2[108];
tmpQN1[187] = + tmpQN2[109];
tmpQN1[188] = + tmpQN2[110];
tmpQN1[189] = + tmpQN2[111];
tmpQN1[190] = 0.0;
;
tmpQN1[191] = 0.0;
;
tmpQN1[192] = 0.0;
;
tmpQN1[193] = 0.0;
;
tmpQN1[194] = 0.0;
;
tmpQN1[195] = 0.0;
;
}

void nmpc_evaluateObjective(  )
{
int runObj;
for (runObj = 0; runObj < 50; ++runObj)
{
nmpcWorkspace.objValueIn[0] = nmpcVariables.x[runObj * 14];
nmpcWorkspace.objValueIn[1] = nmpcVariables.x[runObj * 14 + 1];
nmpcWorkspace.objValueIn[2] = nmpcVariables.x[runObj * 14 + 2];
nmpcWorkspace.objValueIn[3] = nmpcVariables.x[runObj * 14 + 3];
nmpcWorkspace.objValueIn[4] = nmpcVariables.x[runObj * 14 + 4];
nmpcWorkspace.objValueIn[5] = nmpcVariables.x[runObj * 14 + 5];
nmpcWorkspace.objValueIn[6] = nmpcVariables.x[runObj * 14 + 6];
nmpcWorkspace.objValueIn[7] = nmpcVariables.x[runObj * 14 + 7];
nmpcWorkspace.objValueIn[8] = nmpcVariables.x[runObj * 14 + 8];
nmpcWorkspace.objValueIn[9] = nmpcVariables.x[runObj * 14 + 9];
nmpcWorkspace.objValueIn[10] = nmpcVariables.x[runObj * 14 + 10];
nmpcWorkspace.objValueIn[11] = nmpcVariables.x[runObj * 14 + 11];
nmpcWorkspace.objValueIn[12] = nmpcVariables.x[runObj * 14 + 12];
nmpcWorkspace.objValueIn[13] = nmpcVariables.x[runObj * 14 + 13];
nmpcWorkspace.objValueIn[14] = nmpcVariables.u[runObj * 4];
nmpcWorkspace.objValueIn[15] = nmpcVariables.u[runObj * 4 + 1];
nmpcWorkspace.objValueIn[16] = nmpcVariables.u[runObj * 4 + 2];
nmpcWorkspace.objValueIn[17] = nmpcVariables.u[runObj * 4 + 3];
nmpcWorkspace.objValueIn[18] = nmpcVariables.od[runObj * 6];
nmpcWorkspace.objValueIn[19] = nmpcVariables.od[runObj * 6 + 1];
nmpcWorkspace.objValueIn[20] = nmpcVariables.od[runObj * 6 + 2];
nmpcWorkspace.objValueIn[21] = nmpcVariables.od[runObj * 6 + 3];
nmpcWorkspace.objValueIn[22] = nmpcVariables.od[runObj * 6 + 4];
nmpcWorkspace.objValueIn[23] = nmpcVariables.od[runObj * 6 + 5];

nmpc_evaluateLSQ( nmpcWorkspace.objValueIn, nmpcWorkspace.objValueOut );
nmpcWorkspace.Dy[runObj * 13] = nmpcWorkspace.objValueOut[0];
nmpcWorkspace.Dy[runObj * 13 + 1] = nmpcWorkspace.objValueOut[1];
nmpcWorkspace.Dy[runObj * 13 + 2] = nmpcWorkspace.objValueOut[2];
nmpcWorkspace.Dy[runObj * 13 + 3] = nmpcWorkspace.objValueOut[3];
nmpcWorkspace.Dy[runObj * 13 + 4] = nmpcWorkspace.objValueOut[4];
nmpcWorkspace.Dy[runObj * 13 + 5] = nmpcWorkspace.objValueOut[5];
nmpcWorkspace.Dy[runObj * 13 + 6] = nmpcWorkspace.objValueOut[6];
nmpcWorkspace.Dy[runObj * 13 + 7] = nmpcWorkspace.objValueOut[7];
nmpcWorkspace.Dy[runObj * 13 + 8] = nmpcWorkspace.objValueOut[8];
nmpcWorkspace.Dy[runObj * 13 + 9] = nmpcWorkspace.objValueOut[9];
nmpcWorkspace.Dy[runObj * 13 + 10] = nmpcWorkspace.objValueOut[10];
nmpcWorkspace.Dy[runObj * 13 + 11] = nmpcWorkspace.objValueOut[11];
nmpcWorkspace.Dy[runObj * 13 + 12] = nmpcWorkspace.objValueOut[12];

nmpc_setObjQ1Q2( &(nmpcWorkspace.objValueOut[ 13 ]), nmpcVariables.W, &(nmpcWorkspace.Q1[ runObj * 196 ]), &(nmpcWorkspace.Q2[ runObj * 182 ]) );

nmpc_setObjR1R2( nmpcVariables.W, &(nmpcWorkspace.R1[ runObj * 16 ]), &(nmpcWorkspace.R2[ runObj * 52 ]) );

}
nmpcWorkspace.objValueIn[0] = nmpcVariables.x[700];
nmpcWorkspace.objValueIn[1] = nmpcVariables.x[701];
nmpcWorkspace.objValueIn[2] = nmpcVariables.x[702];
nmpcWorkspace.objValueIn[3] = nmpcVariables.x[703];
nmpcWorkspace.objValueIn[4] = nmpcVariables.x[704];
nmpcWorkspace.objValueIn[5] = nmpcVariables.x[705];
nmpcWorkspace.objValueIn[6] = nmpcVariables.x[706];
nmpcWorkspace.objValueIn[7] = nmpcVariables.x[707];
nmpcWorkspace.objValueIn[8] = nmpcVariables.x[708];
nmpcWorkspace.objValueIn[9] = nmpcVariables.x[709];
nmpcWorkspace.objValueIn[10] = nmpcVariables.x[710];
nmpcWorkspace.objValueIn[11] = nmpcVariables.x[711];
nmpcWorkspace.objValueIn[12] = nmpcVariables.x[712];
nmpcWorkspace.objValueIn[13] = nmpcVariables.x[713];
nmpcWorkspace.objValueIn[14] = nmpcVariables.od[300];
nmpcWorkspace.objValueIn[15] = nmpcVariables.od[301];
nmpcWorkspace.objValueIn[16] = nmpcVariables.od[302];
nmpcWorkspace.objValueIn[17] = nmpcVariables.od[303];
nmpcWorkspace.objValueIn[18] = nmpcVariables.od[304];
nmpcWorkspace.objValueIn[19] = nmpcVariables.od[305];
nmpc_evaluateLSQEndTerm( nmpcWorkspace.objValueIn, nmpcWorkspace.objValueOut );

nmpcWorkspace.DyN[0] = nmpcWorkspace.objValueOut[0];
nmpcWorkspace.DyN[1] = nmpcWorkspace.objValueOut[1];
nmpcWorkspace.DyN[2] = nmpcWorkspace.objValueOut[2];
nmpcWorkspace.DyN[3] = nmpcWorkspace.objValueOut[3];
nmpcWorkspace.DyN[4] = nmpcWorkspace.objValueOut[4];
nmpcWorkspace.DyN[5] = nmpcWorkspace.objValueOut[5];
nmpcWorkspace.DyN[6] = nmpcWorkspace.objValueOut[6];
nmpcWorkspace.DyN[7] = nmpcWorkspace.objValueOut[7];

nmpc_setObjQN1QN2( nmpcVariables.WN, nmpcWorkspace.QN1, nmpcWorkspace.QN2 );

}

void nmpc_multGxd( real_t* const dOld, real_t* const Gx1, real_t* const dNew )
{
dNew[0] += + Gx1[0]*dOld[0] + Gx1[1]*dOld[1] + Gx1[2]*dOld[2] + Gx1[3]*dOld[3] + Gx1[4]*dOld[4] + Gx1[5]*dOld[5] + Gx1[6]*dOld[6] + Gx1[7]*dOld[7] + Gx1[8]*dOld[8] + Gx1[9]*dOld[9] + Gx1[10]*dOld[10] + Gx1[11]*dOld[11] + Gx1[12]*dOld[12] + Gx1[13]*dOld[13];
dNew[1] += + Gx1[14]*dOld[0] + Gx1[15]*dOld[1] + Gx1[16]*dOld[2] + Gx1[17]*dOld[3] + Gx1[18]*dOld[4] + Gx1[19]*dOld[5] + Gx1[20]*dOld[6] + Gx1[21]*dOld[7] + Gx1[22]*dOld[8] + Gx1[23]*dOld[9] + Gx1[24]*dOld[10] + Gx1[25]*dOld[11] + Gx1[26]*dOld[12] + Gx1[27]*dOld[13];
dNew[2] += + Gx1[28]*dOld[0] + Gx1[29]*dOld[1] + Gx1[30]*dOld[2] + Gx1[31]*dOld[3] + Gx1[32]*dOld[4] + Gx1[33]*dOld[5] + Gx1[34]*dOld[6] + Gx1[35]*dOld[7] + Gx1[36]*dOld[8] + Gx1[37]*dOld[9] + Gx1[38]*dOld[10] + Gx1[39]*dOld[11] + Gx1[40]*dOld[12] + Gx1[41]*dOld[13];
dNew[3] += + Gx1[42]*dOld[0] + Gx1[43]*dOld[1] + Gx1[44]*dOld[2] + Gx1[45]*dOld[3] + Gx1[46]*dOld[4] + Gx1[47]*dOld[5] + Gx1[48]*dOld[6] + Gx1[49]*dOld[7] + Gx1[50]*dOld[8] + Gx1[51]*dOld[9] + Gx1[52]*dOld[10] + Gx1[53]*dOld[11] + Gx1[54]*dOld[12] + Gx1[55]*dOld[13];
dNew[4] += + Gx1[56]*dOld[0] + Gx1[57]*dOld[1] + Gx1[58]*dOld[2] + Gx1[59]*dOld[3] + Gx1[60]*dOld[4] + Gx1[61]*dOld[5] + Gx1[62]*dOld[6] + Gx1[63]*dOld[7] + Gx1[64]*dOld[8] + Gx1[65]*dOld[9] + Gx1[66]*dOld[10] + Gx1[67]*dOld[11] + Gx1[68]*dOld[12] + Gx1[69]*dOld[13];
dNew[5] += + Gx1[70]*dOld[0] + Gx1[71]*dOld[1] + Gx1[72]*dOld[2] + Gx1[73]*dOld[3] + Gx1[74]*dOld[4] + Gx1[75]*dOld[5] + Gx1[76]*dOld[6] + Gx1[77]*dOld[7] + Gx1[78]*dOld[8] + Gx1[79]*dOld[9] + Gx1[80]*dOld[10] + Gx1[81]*dOld[11] + Gx1[82]*dOld[12] + Gx1[83]*dOld[13];
dNew[6] += + Gx1[84]*dOld[0] + Gx1[85]*dOld[1] + Gx1[86]*dOld[2] + Gx1[87]*dOld[3] + Gx1[88]*dOld[4] + Gx1[89]*dOld[5] + Gx1[90]*dOld[6] + Gx1[91]*dOld[7] + Gx1[92]*dOld[8] + Gx1[93]*dOld[9] + Gx1[94]*dOld[10] + Gx1[95]*dOld[11] + Gx1[96]*dOld[12] + Gx1[97]*dOld[13];
dNew[7] += + Gx1[98]*dOld[0] + Gx1[99]*dOld[1] + Gx1[100]*dOld[2] + Gx1[101]*dOld[3] + Gx1[102]*dOld[4] + Gx1[103]*dOld[5] + Gx1[104]*dOld[6] + Gx1[105]*dOld[7] + Gx1[106]*dOld[8] + Gx1[107]*dOld[9] + Gx1[108]*dOld[10] + Gx1[109]*dOld[11] + Gx1[110]*dOld[12] + Gx1[111]*dOld[13];
dNew[8] += + Gx1[112]*dOld[0] + Gx1[113]*dOld[1] + Gx1[114]*dOld[2] + Gx1[115]*dOld[3] + Gx1[116]*dOld[4] + Gx1[117]*dOld[5] + Gx1[118]*dOld[6] + Gx1[119]*dOld[7] + Gx1[120]*dOld[8] + Gx1[121]*dOld[9] + Gx1[122]*dOld[10] + Gx1[123]*dOld[11] + Gx1[124]*dOld[12] + Gx1[125]*dOld[13];
dNew[9] += + Gx1[126]*dOld[0] + Gx1[127]*dOld[1] + Gx1[128]*dOld[2] + Gx1[129]*dOld[3] + Gx1[130]*dOld[4] + Gx1[131]*dOld[5] + Gx1[132]*dOld[6] + Gx1[133]*dOld[7] + Gx1[134]*dOld[8] + Gx1[135]*dOld[9] + Gx1[136]*dOld[10] + Gx1[137]*dOld[11] + Gx1[138]*dOld[12] + Gx1[139]*dOld[13];
dNew[10] += + Gx1[140]*dOld[0] + Gx1[141]*dOld[1] + Gx1[142]*dOld[2] + Gx1[143]*dOld[3] + Gx1[144]*dOld[4] + Gx1[145]*dOld[5] + Gx1[146]*dOld[6] + Gx1[147]*dOld[7] + Gx1[148]*dOld[8] + Gx1[149]*dOld[9] + Gx1[150]*dOld[10] + Gx1[151]*dOld[11] + Gx1[152]*dOld[12] + Gx1[153]*dOld[13];
dNew[11] += + Gx1[154]*dOld[0] + Gx1[155]*dOld[1] + Gx1[156]*dOld[2] + Gx1[157]*dOld[3] + Gx1[158]*dOld[4] + Gx1[159]*dOld[5] + Gx1[160]*dOld[6] + Gx1[161]*dOld[7] + Gx1[162]*dOld[8] + Gx1[163]*dOld[9] + Gx1[164]*dOld[10] + Gx1[165]*dOld[11] + Gx1[166]*dOld[12] + Gx1[167]*dOld[13];
dNew[12] += + Gx1[168]*dOld[0] + Gx1[169]*dOld[1] + Gx1[170]*dOld[2] + Gx1[171]*dOld[3] + Gx1[172]*dOld[4] + Gx1[173]*dOld[5] + Gx1[174]*dOld[6] + Gx1[175]*dOld[7] + Gx1[176]*dOld[8] + Gx1[177]*dOld[9] + Gx1[178]*dOld[10] + Gx1[179]*dOld[11] + Gx1[180]*dOld[12] + Gx1[181]*dOld[13];
dNew[13] += + Gx1[182]*dOld[0] + Gx1[183]*dOld[1] + Gx1[184]*dOld[2] + Gx1[185]*dOld[3] + Gx1[186]*dOld[4] + Gx1[187]*dOld[5] + Gx1[188]*dOld[6] + Gx1[189]*dOld[7] + Gx1[190]*dOld[8] + Gx1[191]*dOld[9] + Gx1[192]*dOld[10] + Gx1[193]*dOld[11] + Gx1[194]*dOld[12] + Gx1[195]*dOld[13];
}

void nmpc_moveGxT( real_t* const Gx1, real_t* const Gx2 )
{
int lRun1;
int lRun2;
for (lRun1 = 0;lRun1 < 14; ++lRun1)
for (lRun2 = 0;lRun2 < 14; ++lRun2)
Gx2[(lRun1 * 14) + (lRun2)] = Gx1[(lRun1 * 14) + (lRun2)];
}

void nmpc_multGxGx( real_t* const Gx1, real_t* const Gx2, real_t* const Gx3 )
{
Gx3[0] = + Gx1[0]*Gx2[0] + Gx1[1]*Gx2[14] + Gx1[2]*Gx2[28] + Gx1[3]*Gx2[42] + Gx1[4]*Gx2[56] + Gx1[5]*Gx2[70] + Gx1[6]*Gx2[84] + Gx1[7]*Gx2[98] + Gx1[8]*Gx2[112] + Gx1[9]*Gx2[126] + Gx1[10]*Gx2[140] + Gx1[11]*Gx2[154] + Gx1[12]*Gx2[168] + Gx1[13]*Gx2[182];
Gx3[1] = + Gx1[0]*Gx2[1] + Gx1[1]*Gx2[15] + Gx1[2]*Gx2[29] + Gx1[3]*Gx2[43] + Gx1[4]*Gx2[57] + Gx1[5]*Gx2[71] + Gx1[6]*Gx2[85] + Gx1[7]*Gx2[99] + Gx1[8]*Gx2[113] + Gx1[9]*Gx2[127] + Gx1[10]*Gx2[141] + Gx1[11]*Gx2[155] + Gx1[12]*Gx2[169] + Gx1[13]*Gx2[183];
Gx3[2] = + Gx1[0]*Gx2[2] + Gx1[1]*Gx2[16] + Gx1[2]*Gx2[30] + Gx1[3]*Gx2[44] + Gx1[4]*Gx2[58] + Gx1[5]*Gx2[72] + Gx1[6]*Gx2[86] + Gx1[7]*Gx2[100] + Gx1[8]*Gx2[114] + Gx1[9]*Gx2[128] + Gx1[10]*Gx2[142] + Gx1[11]*Gx2[156] + Gx1[12]*Gx2[170] + Gx1[13]*Gx2[184];
Gx3[3] = + Gx1[0]*Gx2[3] + Gx1[1]*Gx2[17] + Gx1[2]*Gx2[31] + Gx1[3]*Gx2[45] + Gx1[4]*Gx2[59] + Gx1[5]*Gx2[73] + Gx1[6]*Gx2[87] + Gx1[7]*Gx2[101] + Gx1[8]*Gx2[115] + Gx1[9]*Gx2[129] + Gx1[10]*Gx2[143] + Gx1[11]*Gx2[157] + Gx1[12]*Gx2[171] + Gx1[13]*Gx2[185];
Gx3[4] = + Gx1[0]*Gx2[4] + Gx1[1]*Gx2[18] + Gx1[2]*Gx2[32] + Gx1[3]*Gx2[46] + Gx1[4]*Gx2[60] + Gx1[5]*Gx2[74] + Gx1[6]*Gx2[88] + Gx1[7]*Gx2[102] + Gx1[8]*Gx2[116] + Gx1[9]*Gx2[130] + Gx1[10]*Gx2[144] + Gx1[11]*Gx2[158] + Gx1[12]*Gx2[172] + Gx1[13]*Gx2[186];
Gx3[5] = + Gx1[0]*Gx2[5] + Gx1[1]*Gx2[19] + Gx1[2]*Gx2[33] + Gx1[3]*Gx2[47] + Gx1[4]*Gx2[61] + Gx1[5]*Gx2[75] + Gx1[6]*Gx2[89] + Gx1[7]*Gx2[103] + Gx1[8]*Gx2[117] + Gx1[9]*Gx2[131] + Gx1[10]*Gx2[145] + Gx1[11]*Gx2[159] + Gx1[12]*Gx2[173] + Gx1[13]*Gx2[187];
Gx3[6] = + Gx1[0]*Gx2[6] + Gx1[1]*Gx2[20] + Gx1[2]*Gx2[34] + Gx1[3]*Gx2[48] + Gx1[4]*Gx2[62] + Gx1[5]*Gx2[76] + Gx1[6]*Gx2[90] + Gx1[7]*Gx2[104] + Gx1[8]*Gx2[118] + Gx1[9]*Gx2[132] + Gx1[10]*Gx2[146] + Gx1[11]*Gx2[160] + Gx1[12]*Gx2[174] + Gx1[13]*Gx2[188];
Gx3[7] = + Gx1[0]*Gx2[7] + Gx1[1]*Gx2[21] + Gx1[2]*Gx2[35] + Gx1[3]*Gx2[49] + Gx1[4]*Gx2[63] + Gx1[5]*Gx2[77] + Gx1[6]*Gx2[91] + Gx1[7]*Gx2[105] + Gx1[8]*Gx2[119] + Gx1[9]*Gx2[133] + Gx1[10]*Gx2[147] + Gx1[11]*Gx2[161] + Gx1[12]*Gx2[175] + Gx1[13]*Gx2[189];
Gx3[8] = + Gx1[0]*Gx2[8] + Gx1[1]*Gx2[22] + Gx1[2]*Gx2[36] + Gx1[3]*Gx2[50] + Gx1[4]*Gx2[64] + Gx1[5]*Gx2[78] + Gx1[6]*Gx2[92] + Gx1[7]*Gx2[106] + Gx1[8]*Gx2[120] + Gx1[9]*Gx2[134] + Gx1[10]*Gx2[148] + Gx1[11]*Gx2[162] + Gx1[12]*Gx2[176] + Gx1[13]*Gx2[190];
Gx3[9] = + Gx1[0]*Gx2[9] + Gx1[1]*Gx2[23] + Gx1[2]*Gx2[37] + Gx1[3]*Gx2[51] + Gx1[4]*Gx2[65] + Gx1[5]*Gx2[79] + Gx1[6]*Gx2[93] + Gx1[7]*Gx2[107] + Gx1[8]*Gx2[121] + Gx1[9]*Gx2[135] + Gx1[10]*Gx2[149] + Gx1[11]*Gx2[163] + Gx1[12]*Gx2[177] + Gx1[13]*Gx2[191];
Gx3[10] = + Gx1[0]*Gx2[10] + Gx1[1]*Gx2[24] + Gx1[2]*Gx2[38] + Gx1[3]*Gx2[52] + Gx1[4]*Gx2[66] + Gx1[5]*Gx2[80] + Gx1[6]*Gx2[94] + Gx1[7]*Gx2[108] + Gx1[8]*Gx2[122] + Gx1[9]*Gx2[136] + Gx1[10]*Gx2[150] + Gx1[11]*Gx2[164] + Gx1[12]*Gx2[178] + Gx1[13]*Gx2[192];
Gx3[11] = + Gx1[0]*Gx2[11] + Gx1[1]*Gx2[25] + Gx1[2]*Gx2[39] + Gx1[3]*Gx2[53] + Gx1[4]*Gx2[67] + Gx1[5]*Gx2[81] + Gx1[6]*Gx2[95] + Gx1[7]*Gx2[109] + Gx1[8]*Gx2[123] + Gx1[9]*Gx2[137] + Gx1[10]*Gx2[151] + Gx1[11]*Gx2[165] + Gx1[12]*Gx2[179] + Gx1[13]*Gx2[193];
Gx3[12] = + Gx1[0]*Gx2[12] + Gx1[1]*Gx2[26] + Gx1[2]*Gx2[40] + Gx1[3]*Gx2[54] + Gx1[4]*Gx2[68] + Gx1[5]*Gx2[82] + Gx1[6]*Gx2[96] + Gx1[7]*Gx2[110] + Gx1[8]*Gx2[124] + Gx1[9]*Gx2[138] + Gx1[10]*Gx2[152] + Gx1[11]*Gx2[166] + Gx1[12]*Gx2[180] + Gx1[13]*Gx2[194];
Gx3[13] = + Gx1[0]*Gx2[13] + Gx1[1]*Gx2[27] + Gx1[2]*Gx2[41] + Gx1[3]*Gx2[55] + Gx1[4]*Gx2[69] + Gx1[5]*Gx2[83] + Gx1[6]*Gx2[97] + Gx1[7]*Gx2[111] + Gx1[8]*Gx2[125] + Gx1[9]*Gx2[139] + Gx1[10]*Gx2[153] + Gx1[11]*Gx2[167] + Gx1[12]*Gx2[181] + Gx1[13]*Gx2[195];
Gx3[14] = + Gx1[14]*Gx2[0] + Gx1[15]*Gx2[14] + Gx1[16]*Gx2[28] + Gx1[17]*Gx2[42] + Gx1[18]*Gx2[56] + Gx1[19]*Gx2[70] + Gx1[20]*Gx2[84] + Gx1[21]*Gx2[98] + Gx1[22]*Gx2[112] + Gx1[23]*Gx2[126] + Gx1[24]*Gx2[140] + Gx1[25]*Gx2[154] + Gx1[26]*Gx2[168] + Gx1[27]*Gx2[182];
Gx3[15] = + Gx1[14]*Gx2[1] + Gx1[15]*Gx2[15] + Gx1[16]*Gx2[29] + Gx1[17]*Gx2[43] + Gx1[18]*Gx2[57] + Gx1[19]*Gx2[71] + Gx1[20]*Gx2[85] + Gx1[21]*Gx2[99] + Gx1[22]*Gx2[113] + Gx1[23]*Gx2[127] + Gx1[24]*Gx2[141] + Gx1[25]*Gx2[155] + Gx1[26]*Gx2[169] + Gx1[27]*Gx2[183];
Gx3[16] = + Gx1[14]*Gx2[2] + Gx1[15]*Gx2[16] + Gx1[16]*Gx2[30] + Gx1[17]*Gx2[44] + Gx1[18]*Gx2[58] + Gx1[19]*Gx2[72] + Gx1[20]*Gx2[86] + Gx1[21]*Gx2[100] + Gx1[22]*Gx2[114] + Gx1[23]*Gx2[128] + Gx1[24]*Gx2[142] + Gx1[25]*Gx2[156] + Gx1[26]*Gx2[170] + Gx1[27]*Gx2[184];
Gx3[17] = + Gx1[14]*Gx2[3] + Gx1[15]*Gx2[17] + Gx1[16]*Gx2[31] + Gx1[17]*Gx2[45] + Gx1[18]*Gx2[59] + Gx1[19]*Gx2[73] + Gx1[20]*Gx2[87] + Gx1[21]*Gx2[101] + Gx1[22]*Gx2[115] + Gx1[23]*Gx2[129] + Gx1[24]*Gx2[143] + Gx1[25]*Gx2[157] + Gx1[26]*Gx2[171] + Gx1[27]*Gx2[185];
Gx3[18] = + Gx1[14]*Gx2[4] + Gx1[15]*Gx2[18] + Gx1[16]*Gx2[32] + Gx1[17]*Gx2[46] + Gx1[18]*Gx2[60] + Gx1[19]*Gx2[74] + Gx1[20]*Gx2[88] + Gx1[21]*Gx2[102] + Gx1[22]*Gx2[116] + Gx1[23]*Gx2[130] + Gx1[24]*Gx2[144] + Gx1[25]*Gx2[158] + Gx1[26]*Gx2[172] + Gx1[27]*Gx2[186];
Gx3[19] = + Gx1[14]*Gx2[5] + Gx1[15]*Gx2[19] + Gx1[16]*Gx2[33] + Gx1[17]*Gx2[47] + Gx1[18]*Gx2[61] + Gx1[19]*Gx2[75] + Gx1[20]*Gx2[89] + Gx1[21]*Gx2[103] + Gx1[22]*Gx2[117] + Gx1[23]*Gx2[131] + Gx1[24]*Gx2[145] + Gx1[25]*Gx2[159] + Gx1[26]*Gx2[173] + Gx1[27]*Gx2[187];
Gx3[20] = + Gx1[14]*Gx2[6] + Gx1[15]*Gx2[20] + Gx1[16]*Gx2[34] + Gx1[17]*Gx2[48] + Gx1[18]*Gx2[62] + Gx1[19]*Gx2[76] + Gx1[20]*Gx2[90] + Gx1[21]*Gx2[104] + Gx1[22]*Gx2[118] + Gx1[23]*Gx2[132] + Gx1[24]*Gx2[146] + Gx1[25]*Gx2[160] + Gx1[26]*Gx2[174] + Gx1[27]*Gx2[188];
Gx3[21] = + Gx1[14]*Gx2[7] + Gx1[15]*Gx2[21] + Gx1[16]*Gx2[35] + Gx1[17]*Gx2[49] + Gx1[18]*Gx2[63] + Gx1[19]*Gx2[77] + Gx1[20]*Gx2[91] + Gx1[21]*Gx2[105] + Gx1[22]*Gx2[119] + Gx1[23]*Gx2[133] + Gx1[24]*Gx2[147] + Gx1[25]*Gx2[161] + Gx1[26]*Gx2[175] + Gx1[27]*Gx2[189];
Gx3[22] = + Gx1[14]*Gx2[8] + Gx1[15]*Gx2[22] + Gx1[16]*Gx2[36] + Gx1[17]*Gx2[50] + Gx1[18]*Gx2[64] + Gx1[19]*Gx2[78] + Gx1[20]*Gx2[92] + Gx1[21]*Gx2[106] + Gx1[22]*Gx2[120] + Gx1[23]*Gx2[134] + Gx1[24]*Gx2[148] + Gx1[25]*Gx2[162] + Gx1[26]*Gx2[176] + Gx1[27]*Gx2[190];
Gx3[23] = + Gx1[14]*Gx2[9] + Gx1[15]*Gx2[23] + Gx1[16]*Gx2[37] + Gx1[17]*Gx2[51] + Gx1[18]*Gx2[65] + Gx1[19]*Gx2[79] + Gx1[20]*Gx2[93] + Gx1[21]*Gx2[107] + Gx1[22]*Gx2[121] + Gx1[23]*Gx2[135] + Gx1[24]*Gx2[149] + Gx1[25]*Gx2[163] + Gx1[26]*Gx2[177] + Gx1[27]*Gx2[191];
Gx3[24] = + Gx1[14]*Gx2[10] + Gx1[15]*Gx2[24] + Gx1[16]*Gx2[38] + Gx1[17]*Gx2[52] + Gx1[18]*Gx2[66] + Gx1[19]*Gx2[80] + Gx1[20]*Gx2[94] + Gx1[21]*Gx2[108] + Gx1[22]*Gx2[122] + Gx1[23]*Gx2[136] + Gx1[24]*Gx2[150] + Gx1[25]*Gx2[164] + Gx1[26]*Gx2[178] + Gx1[27]*Gx2[192];
Gx3[25] = + Gx1[14]*Gx2[11] + Gx1[15]*Gx2[25] + Gx1[16]*Gx2[39] + Gx1[17]*Gx2[53] + Gx1[18]*Gx2[67] + Gx1[19]*Gx2[81] + Gx1[20]*Gx2[95] + Gx1[21]*Gx2[109] + Gx1[22]*Gx2[123] + Gx1[23]*Gx2[137] + Gx1[24]*Gx2[151] + Gx1[25]*Gx2[165] + Gx1[26]*Gx2[179] + Gx1[27]*Gx2[193];
Gx3[26] = + Gx1[14]*Gx2[12] + Gx1[15]*Gx2[26] + Gx1[16]*Gx2[40] + Gx1[17]*Gx2[54] + Gx1[18]*Gx2[68] + Gx1[19]*Gx2[82] + Gx1[20]*Gx2[96] + Gx1[21]*Gx2[110] + Gx1[22]*Gx2[124] + Gx1[23]*Gx2[138] + Gx1[24]*Gx2[152] + Gx1[25]*Gx2[166] + Gx1[26]*Gx2[180] + Gx1[27]*Gx2[194];
Gx3[27] = + Gx1[14]*Gx2[13] + Gx1[15]*Gx2[27] + Gx1[16]*Gx2[41] + Gx1[17]*Gx2[55] + Gx1[18]*Gx2[69] + Gx1[19]*Gx2[83] + Gx1[20]*Gx2[97] + Gx1[21]*Gx2[111] + Gx1[22]*Gx2[125] + Gx1[23]*Gx2[139] + Gx1[24]*Gx2[153] + Gx1[25]*Gx2[167] + Gx1[26]*Gx2[181] + Gx1[27]*Gx2[195];
Gx3[28] = + Gx1[28]*Gx2[0] + Gx1[29]*Gx2[14] + Gx1[30]*Gx2[28] + Gx1[31]*Gx2[42] + Gx1[32]*Gx2[56] + Gx1[33]*Gx2[70] + Gx1[34]*Gx2[84] + Gx1[35]*Gx2[98] + Gx1[36]*Gx2[112] + Gx1[37]*Gx2[126] + Gx1[38]*Gx2[140] + Gx1[39]*Gx2[154] + Gx1[40]*Gx2[168] + Gx1[41]*Gx2[182];
Gx3[29] = + Gx1[28]*Gx2[1] + Gx1[29]*Gx2[15] + Gx1[30]*Gx2[29] + Gx1[31]*Gx2[43] + Gx1[32]*Gx2[57] + Gx1[33]*Gx2[71] + Gx1[34]*Gx2[85] + Gx1[35]*Gx2[99] + Gx1[36]*Gx2[113] + Gx1[37]*Gx2[127] + Gx1[38]*Gx2[141] + Gx1[39]*Gx2[155] + Gx1[40]*Gx2[169] + Gx1[41]*Gx2[183];
Gx3[30] = + Gx1[28]*Gx2[2] + Gx1[29]*Gx2[16] + Gx1[30]*Gx2[30] + Gx1[31]*Gx2[44] + Gx1[32]*Gx2[58] + Gx1[33]*Gx2[72] + Gx1[34]*Gx2[86] + Gx1[35]*Gx2[100] + Gx1[36]*Gx2[114] + Gx1[37]*Gx2[128] + Gx1[38]*Gx2[142] + Gx1[39]*Gx2[156] + Gx1[40]*Gx2[170] + Gx1[41]*Gx2[184];
Gx3[31] = + Gx1[28]*Gx2[3] + Gx1[29]*Gx2[17] + Gx1[30]*Gx2[31] + Gx1[31]*Gx2[45] + Gx1[32]*Gx2[59] + Gx1[33]*Gx2[73] + Gx1[34]*Gx2[87] + Gx1[35]*Gx2[101] + Gx1[36]*Gx2[115] + Gx1[37]*Gx2[129] + Gx1[38]*Gx2[143] + Gx1[39]*Gx2[157] + Gx1[40]*Gx2[171] + Gx1[41]*Gx2[185];
Gx3[32] = + Gx1[28]*Gx2[4] + Gx1[29]*Gx2[18] + Gx1[30]*Gx2[32] + Gx1[31]*Gx2[46] + Gx1[32]*Gx2[60] + Gx1[33]*Gx2[74] + Gx1[34]*Gx2[88] + Gx1[35]*Gx2[102] + Gx1[36]*Gx2[116] + Gx1[37]*Gx2[130] + Gx1[38]*Gx2[144] + Gx1[39]*Gx2[158] + Gx1[40]*Gx2[172] + Gx1[41]*Gx2[186];
Gx3[33] = + Gx1[28]*Gx2[5] + Gx1[29]*Gx2[19] + Gx1[30]*Gx2[33] + Gx1[31]*Gx2[47] + Gx1[32]*Gx2[61] + Gx1[33]*Gx2[75] + Gx1[34]*Gx2[89] + Gx1[35]*Gx2[103] + Gx1[36]*Gx2[117] + Gx1[37]*Gx2[131] + Gx1[38]*Gx2[145] + Gx1[39]*Gx2[159] + Gx1[40]*Gx2[173] + Gx1[41]*Gx2[187];
Gx3[34] = + Gx1[28]*Gx2[6] + Gx1[29]*Gx2[20] + Gx1[30]*Gx2[34] + Gx1[31]*Gx2[48] + Gx1[32]*Gx2[62] + Gx1[33]*Gx2[76] + Gx1[34]*Gx2[90] + Gx1[35]*Gx2[104] + Gx1[36]*Gx2[118] + Gx1[37]*Gx2[132] + Gx1[38]*Gx2[146] + Gx1[39]*Gx2[160] + Gx1[40]*Gx2[174] + Gx1[41]*Gx2[188];
Gx3[35] = + Gx1[28]*Gx2[7] + Gx1[29]*Gx2[21] + Gx1[30]*Gx2[35] + Gx1[31]*Gx2[49] + Gx1[32]*Gx2[63] + Gx1[33]*Gx2[77] + Gx1[34]*Gx2[91] + Gx1[35]*Gx2[105] + Gx1[36]*Gx2[119] + Gx1[37]*Gx2[133] + Gx1[38]*Gx2[147] + Gx1[39]*Gx2[161] + Gx1[40]*Gx2[175] + Gx1[41]*Gx2[189];
Gx3[36] = + Gx1[28]*Gx2[8] + Gx1[29]*Gx2[22] + Gx1[30]*Gx2[36] + Gx1[31]*Gx2[50] + Gx1[32]*Gx2[64] + Gx1[33]*Gx2[78] + Gx1[34]*Gx2[92] + Gx1[35]*Gx2[106] + Gx1[36]*Gx2[120] + Gx1[37]*Gx2[134] + Gx1[38]*Gx2[148] + Gx1[39]*Gx2[162] + Gx1[40]*Gx2[176] + Gx1[41]*Gx2[190];
Gx3[37] = + Gx1[28]*Gx2[9] + Gx1[29]*Gx2[23] + Gx1[30]*Gx2[37] + Gx1[31]*Gx2[51] + Gx1[32]*Gx2[65] + Gx1[33]*Gx2[79] + Gx1[34]*Gx2[93] + Gx1[35]*Gx2[107] + Gx1[36]*Gx2[121] + Gx1[37]*Gx2[135] + Gx1[38]*Gx2[149] + Gx1[39]*Gx2[163] + Gx1[40]*Gx2[177] + Gx1[41]*Gx2[191];
Gx3[38] = + Gx1[28]*Gx2[10] + Gx1[29]*Gx2[24] + Gx1[30]*Gx2[38] + Gx1[31]*Gx2[52] + Gx1[32]*Gx2[66] + Gx1[33]*Gx2[80] + Gx1[34]*Gx2[94] + Gx1[35]*Gx2[108] + Gx1[36]*Gx2[122] + Gx1[37]*Gx2[136] + Gx1[38]*Gx2[150] + Gx1[39]*Gx2[164] + Gx1[40]*Gx2[178] + Gx1[41]*Gx2[192];
Gx3[39] = + Gx1[28]*Gx2[11] + Gx1[29]*Gx2[25] + Gx1[30]*Gx2[39] + Gx1[31]*Gx2[53] + Gx1[32]*Gx2[67] + Gx1[33]*Gx2[81] + Gx1[34]*Gx2[95] + Gx1[35]*Gx2[109] + Gx1[36]*Gx2[123] + Gx1[37]*Gx2[137] + Gx1[38]*Gx2[151] + Gx1[39]*Gx2[165] + Gx1[40]*Gx2[179] + Gx1[41]*Gx2[193];
Gx3[40] = + Gx1[28]*Gx2[12] + Gx1[29]*Gx2[26] + Gx1[30]*Gx2[40] + Gx1[31]*Gx2[54] + Gx1[32]*Gx2[68] + Gx1[33]*Gx2[82] + Gx1[34]*Gx2[96] + Gx1[35]*Gx2[110] + Gx1[36]*Gx2[124] + Gx1[37]*Gx2[138] + Gx1[38]*Gx2[152] + Gx1[39]*Gx2[166] + Gx1[40]*Gx2[180] + Gx1[41]*Gx2[194];
Gx3[41] = + Gx1[28]*Gx2[13] + Gx1[29]*Gx2[27] + Gx1[30]*Gx2[41] + Gx1[31]*Gx2[55] + Gx1[32]*Gx2[69] + Gx1[33]*Gx2[83] + Gx1[34]*Gx2[97] + Gx1[35]*Gx2[111] + Gx1[36]*Gx2[125] + Gx1[37]*Gx2[139] + Gx1[38]*Gx2[153] + Gx1[39]*Gx2[167] + Gx1[40]*Gx2[181] + Gx1[41]*Gx2[195];
Gx3[42] = + Gx1[42]*Gx2[0] + Gx1[43]*Gx2[14] + Gx1[44]*Gx2[28] + Gx1[45]*Gx2[42] + Gx1[46]*Gx2[56] + Gx1[47]*Gx2[70] + Gx1[48]*Gx2[84] + Gx1[49]*Gx2[98] + Gx1[50]*Gx2[112] + Gx1[51]*Gx2[126] + Gx1[52]*Gx2[140] + Gx1[53]*Gx2[154] + Gx1[54]*Gx2[168] + Gx1[55]*Gx2[182];
Gx3[43] = + Gx1[42]*Gx2[1] + Gx1[43]*Gx2[15] + Gx1[44]*Gx2[29] + Gx1[45]*Gx2[43] + Gx1[46]*Gx2[57] + Gx1[47]*Gx2[71] + Gx1[48]*Gx2[85] + Gx1[49]*Gx2[99] + Gx1[50]*Gx2[113] + Gx1[51]*Gx2[127] + Gx1[52]*Gx2[141] + Gx1[53]*Gx2[155] + Gx1[54]*Gx2[169] + Gx1[55]*Gx2[183];
Gx3[44] = + Gx1[42]*Gx2[2] + Gx1[43]*Gx2[16] + Gx1[44]*Gx2[30] + Gx1[45]*Gx2[44] + Gx1[46]*Gx2[58] + Gx1[47]*Gx2[72] + Gx1[48]*Gx2[86] + Gx1[49]*Gx2[100] + Gx1[50]*Gx2[114] + Gx1[51]*Gx2[128] + Gx1[52]*Gx2[142] + Gx1[53]*Gx2[156] + Gx1[54]*Gx2[170] + Gx1[55]*Gx2[184];
Gx3[45] = + Gx1[42]*Gx2[3] + Gx1[43]*Gx2[17] + Gx1[44]*Gx2[31] + Gx1[45]*Gx2[45] + Gx1[46]*Gx2[59] + Gx1[47]*Gx2[73] + Gx1[48]*Gx2[87] + Gx1[49]*Gx2[101] + Gx1[50]*Gx2[115] + Gx1[51]*Gx2[129] + Gx1[52]*Gx2[143] + Gx1[53]*Gx2[157] + Gx1[54]*Gx2[171] + Gx1[55]*Gx2[185];
Gx3[46] = + Gx1[42]*Gx2[4] + Gx1[43]*Gx2[18] + Gx1[44]*Gx2[32] + Gx1[45]*Gx2[46] + Gx1[46]*Gx2[60] + Gx1[47]*Gx2[74] + Gx1[48]*Gx2[88] + Gx1[49]*Gx2[102] + Gx1[50]*Gx2[116] + Gx1[51]*Gx2[130] + Gx1[52]*Gx2[144] + Gx1[53]*Gx2[158] + Gx1[54]*Gx2[172] + Gx1[55]*Gx2[186];
Gx3[47] = + Gx1[42]*Gx2[5] + Gx1[43]*Gx2[19] + Gx1[44]*Gx2[33] + Gx1[45]*Gx2[47] + Gx1[46]*Gx2[61] + Gx1[47]*Gx2[75] + Gx1[48]*Gx2[89] + Gx1[49]*Gx2[103] + Gx1[50]*Gx2[117] + Gx1[51]*Gx2[131] + Gx1[52]*Gx2[145] + Gx1[53]*Gx2[159] + Gx1[54]*Gx2[173] + Gx1[55]*Gx2[187];
Gx3[48] = + Gx1[42]*Gx2[6] + Gx1[43]*Gx2[20] + Gx1[44]*Gx2[34] + Gx1[45]*Gx2[48] + Gx1[46]*Gx2[62] + Gx1[47]*Gx2[76] + Gx1[48]*Gx2[90] + Gx1[49]*Gx2[104] + Gx1[50]*Gx2[118] + Gx1[51]*Gx2[132] + Gx1[52]*Gx2[146] + Gx1[53]*Gx2[160] + Gx1[54]*Gx2[174] + Gx1[55]*Gx2[188];
Gx3[49] = + Gx1[42]*Gx2[7] + Gx1[43]*Gx2[21] + Gx1[44]*Gx2[35] + Gx1[45]*Gx2[49] + Gx1[46]*Gx2[63] + Gx1[47]*Gx2[77] + Gx1[48]*Gx2[91] + Gx1[49]*Gx2[105] + Gx1[50]*Gx2[119] + Gx1[51]*Gx2[133] + Gx1[52]*Gx2[147] + Gx1[53]*Gx2[161] + Gx1[54]*Gx2[175] + Gx1[55]*Gx2[189];
Gx3[50] = + Gx1[42]*Gx2[8] + Gx1[43]*Gx2[22] + Gx1[44]*Gx2[36] + Gx1[45]*Gx2[50] + Gx1[46]*Gx2[64] + Gx1[47]*Gx2[78] + Gx1[48]*Gx2[92] + Gx1[49]*Gx2[106] + Gx1[50]*Gx2[120] + Gx1[51]*Gx2[134] + Gx1[52]*Gx2[148] + Gx1[53]*Gx2[162] + Gx1[54]*Gx2[176] + Gx1[55]*Gx2[190];
Gx3[51] = + Gx1[42]*Gx2[9] + Gx1[43]*Gx2[23] + Gx1[44]*Gx2[37] + Gx1[45]*Gx2[51] + Gx1[46]*Gx2[65] + Gx1[47]*Gx2[79] + Gx1[48]*Gx2[93] + Gx1[49]*Gx2[107] + Gx1[50]*Gx2[121] + Gx1[51]*Gx2[135] + Gx1[52]*Gx2[149] + Gx1[53]*Gx2[163] + Gx1[54]*Gx2[177] + Gx1[55]*Gx2[191];
Gx3[52] = + Gx1[42]*Gx2[10] + Gx1[43]*Gx2[24] + Gx1[44]*Gx2[38] + Gx1[45]*Gx2[52] + Gx1[46]*Gx2[66] + Gx1[47]*Gx2[80] + Gx1[48]*Gx2[94] + Gx1[49]*Gx2[108] + Gx1[50]*Gx2[122] + Gx1[51]*Gx2[136] + Gx1[52]*Gx2[150] + Gx1[53]*Gx2[164] + Gx1[54]*Gx2[178] + Gx1[55]*Gx2[192];
Gx3[53] = + Gx1[42]*Gx2[11] + Gx1[43]*Gx2[25] + Gx1[44]*Gx2[39] + Gx1[45]*Gx2[53] + Gx1[46]*Gx2[67] + Gx1[47]*Gx2[81] + Gx1[48]*Gx2[95] + Gx1[49]*Gx2[109] + Gx1[50]*Gx2[123] + Gx1[51]*Gx2[137] + Gx1[52]*Gx2[151] + Gx1[53]*Gx2[165] + Gx1[54]*Gx2[179] + Gx1[55]*Gx2[193];
Gx3[54] = + Gx1[42]*Gx2[12] + Gx1[43]*Gx2[26] + Gx1[44]*Gx2[40] + Gx1[45]*Gx2[54] + Gx1[46]*Gx2[68] + Gx1[47]*Gx2[82] + Gx1[48]*Gx2[96] + Gx1[49]*Gx2[110] + Gx1[50]*Gx2[124] + Gx1[51]*Gx2[138] + Gx1[52]*Gx2[152] + Gx1[53]*Gx2[166] + Gx1[54]*Gx2[180] + Gx1[55]*Gx2[194];
Gx3[55] = + Gx1[42]*Gx2[13] + Gx1[43]*Gx2[27] + Gx1[44]*Gx2[41] + Gx1[45]*Gx2[55] + Gx1[46]*Gx2[69] + Gx1[47]*Gx2[83] + Gx1[48]*Gx2[97] + Gx1[49]*Gx2[111] + Gx1[50]*Gx2[125] + Gx1[51]*Gx2[139] + Gx1[52]*Gx2[153] + Gx1[53]*Gx2[167] + Gx1[54]*Gx2[181] + Gx1[55]*Gx2[195];
Gx3[56] = + Gx1[56]*Gx2[0] + Gx1[57]*Gx2[14] + Gx1[58]*Gx2[28] + Gx1[59]*Gx2[42] + Gx1[60]*Gx2[56] + Gx1[61]*Gx2[70] + Gx1[62]*Gx2[84] + Gx1[63]*Gx2[98] + Gx1[64]*Gx2[112] + Gx1[65]*Gx2[126] + Gx1[66]*Gx2[140] + Gx1[67]*Gx2[154] + Gx1[68]*Gx2[168] + Gx1[69]*Gx2[182];
Gx3[57] = + Gx1[56]*Gx2[1] + Gx1[57]*Gx2[15] + Gx1[58]*Gx2[29] + Gx1[59]*Gx2[43] + Gx1[60]*Gx2[57] + Gx1[61]*Gx2[71] + Gx1[62]*Gx2[85] + Gx1[63]*Gx2[99] + Gx1[64]*Gx2[113] + Gx1[65]*Gx2[127] + Gx1[66]*Gx2[141] + Gx1[67]*Gx2[155] + Gx1[68]*Gx2[169] + Gx1[69]*Gx2[183];
Gx3[58] = + Gx1[56]*Gx2[2] + Gx1[57]*Gx2[16] + Gx1[58]*Gx2[30] + Gx1[59]*Gx2[44] + Gx1[60]*Gx2[58] + Gx1[61]*Gx2[72] + Gx1[62]*Gx2[86] + Gx1[63]*Gx2[100] + Gx1[64]*Gx2[114] + Gx1[65]*Gx2[128] + Gx1[66]*Gx2[142] + Gx1[67]*Gx2[156] + Gx1[68]*Gx2[170] + Gx1[69]*Gx2[184];
Gx3[59] = + Gx1[56]*Gx2[3] + Gx1[57]*Gx2[17] + Gx1[58]*Gx2[31] + Gx1[59]*Gx2[45] + Gx1[60]*Gx2[59] + Gx1[61]*Gx2[73] + Gx1[62]*Gx2[87] + Gx1[63]*Gx2[101] + Gx1[64]*Gx2[115] + Gx1[65]*Gx2[129] + Gx1[66]*Gx2[143] + Gx1[67]*Gx2[157] + Gx1[68]*Gx2[171] + Gx1[69]*Gx2[185];
Gx3[60] = + Gx1[56]*Gx2[4] + Gx1[57]*Gx2[18] + Gx1[58]*Gx2[32] + Gx1[59]*Gx2[46] + Gx1[60]*Gx2[60] + Gx1[61]*Gx2[74] + Gx1[62]*Gx2[88] + Gx1[63]*Gx2[102] + Gx1[64]*Gx2[116] + Gx1[65]*Gx2[130] + Gx1[66]*Gx2[144] + Gx1[67]*Gx2[158] + Gx1[68]*Gx2[172] + Gx1[69]*Gx2[186];
Gx3[61] = + Gx1[56]*Gx2[5] + Gx1[57]*Gx2[19] + Gx1[58]*Gx2[33] + Gx1[59]*Gx2[47] + Gx1[60]*Gx2[61] + Gx1[61]*Gx2[75] + Gx1[62]*Gx2[89] + Gx1[63]*Gx2[103] + Gx1[64]*Gx2[117] + Gx1[65]*Gx2[131] + Gx1[66]*Gx2[145] + Gx1[67]*Gx2[159] + Gx1[68]*Gx2[173] + Gx1[69]*Gx2[187];
Gx3[62] = + Gx1[56]*Gx2[6] + Gx1[57]*Gx2[20] + Gx1[58]*Gx2[34] + Gx1[59]*Gx2[48] + Gx1[60]*Gx2[62] + Gx1[61]*Gx2[76] + Gx1[62]*Gx2[90] + Gx1[63]*Gx2[104] + Gx1[64]*Gx2[118] + Gx1[65]*Gx2[132] + Gx1[66]*Gx2[146] + Gx1[67]*Gx2[160] + Gx1[68]*Gx2[174] + Gx1[69]*Gx2[188];
Gx3[63] = + Gx1[56]*Gx2[7] + Gx1[57]*Gx2[21] + Gx1[58]*Gx2[35] + Gx1[59]*Gx2[49] + Gx1[60]*Gx2[63] + Gx1[61]*Gx2[77] + Gx1[62]*Gx2[91] + Gx1[63]*Gx2[105] + Gx1[64]*Gx2[119] + Gx1[65]*Gx2[133] + Gx1[66]*Gx2[147] + Gx1[67]*Gx2[161] + Gx1[68]*Gx2[175] + Gx1[69]*Gx2[189];
Gx3[64] = + Gx1[56]*Gx2[8] + Gx1[57]*Gx2[22] + Gx1[58]*Gx2[36] + Gx1[59]*Gx2[50] + Gx1[60]*Gx2[64] + Gx1[61]*Gx2[78] + Gx1[62]*Gx2[92] + Gx1[63]*Gx2[106] + Gx1[64]*Gx2[120] + Gx1[65]*Gx2[134] + Gx1[66]*Gx2[148] + Gx1[67]*Gx2[162] + Gx1[68]*Gx2[176] + Gx1[69]*Gx2[190];
Gx3[65] = + Gx1[56]*Gx2[9] + Gx1[57]*Gx2[23] + Gx1[58]*Gx2[37] + Gx1[59]*Gx2[51] + Gx1[60]*Gx2[65] + Gx1[61]*Gx2[79] + Gx1[62]*Gx2[93] + Gx1[63]*Gx2[107] + Gx1[64]*Gx2[121] + Gx1[65]*Gx2[135] + Gx1[66]*Gx2[149] + Gx1[67]*Gx2[163] + Gx1[68]*Gx2[177] + Gx1[69]*Gx2[191];
Gx3[66] = + Gx1[56]*Gx2[10] + Gx1[57]*Gx2[24] + Gx1[58]*Gx2[38] + Gx1[59]*Gx2[52] + Gx1[60]*Gx2[66] + Gx1[61]*Gx2[80] + Gx1[62]*Gx2[94] + Gx1[63]*Gx2[108] + Gx1[64]*Gx2[122] + Gx1[65]*Gx2[136] + Gx1[66]*Gx2[150] + Gx1[67]*Gx2[164] + Gx1[68]*Gx2[178] + Gx1[69]*Gx2[192];
Gx3[67] = + Gx1[56]*Gx2[11] + Gx1[57]*Gx2[25] + Gx1[58]*Gx2[39] + Gx1[59]*Gx2[53] + Gx1[60]*Gx2[67] + Gx1[61]*Gx2[81] + Gx1[62]*Gx2[95] + Gx1[63]*Gx2[109] + Gx1[64]*Gx2[123] + Gx1[65]*Gx2[137] + Gx1[66]*Gx2[151] + Gx1[67]*Gx2[165] + Gx1[68]*Gx2[179] + Gx1[69]*Gx2[193];
Gx3[68] = + Gx1[56]*Gx2[12] + Gx1[57]*Gx2[26] + Gx1[58]*Gx2[40] + Gx1[59]*Gx2[54] + Gx1[60]*Gx2[68] + Gx1[61]*Gx2[82] + Gx1[62]*Gx2[96] + Gx1[63]*Gx2[110] + Gx1[64]*Gx2[124] + Gx1[65]*Gx2[138] + Gx1[66]*Gx2[152] + Gx1[67]*Gx2[166] + Gx1[68]*Gx2[180] + Gx1[69]*Gx2[194];
Gx3[69] = + Gx1[56]*Gx2[13] + Gx1[57]*Gx2[27] + Gx1[58]*Gx2[41] + Gx1[59]*Gx2[55] + Gx1[60]*Gx2[69] + Gx1[61]*Gx2[83] + Gx1[62]*Gx2[97] + Gx1[63]*Gx2[111] + Gx1[64]*Gx2[125] + Gx1[65]*Gx2[139] + Gx1[66]*Gx2[153] + Gx1[67]*Gx2[167] + Gx1[68]*Gx2[181] + Gx1[69]*Gx2[195];
Gx3[70] = + Gx1[70]*Gx2[0] + Gx1[71]*Gx2[14] + Gx1[72]*Gx2[28] + Gx1[73]*Gx2[42] + Gx1[74]*Gx2[56] + Gx1[75]*Gx2[70] + Gx1[76]*Gx2[84] + Gx1[77]*Gx2[98] + Gx1[78]*Gx2[112] + Gx1[79]*Gx2[126] + Gx1[80]*Gx2[140] + Gx1[81]*Gx2[154] + Gx1[82]*Gx2[168] + Gx1[83]*Gx2[182];
Gx3[71] = + Gx1[70]*Gx2[1] + Gx1[71]*Gx2[15] + Gx1[72]*Gx2[29] + Gx1[73]*Gx2[43] + Gx1[74]*Gx2[57] + Gx1[75]*Gx2[71] + Gx1[76]*Gx2[85] + Gx1[77]*Gx2[99] + Gx1[78]*Gx2[113] + Gx1[79]*Gx2[127] + Gx1[80]*Gx2[141] + Gx1[81]*Gx2[155] + Gx1[82]*Gx2[169] + Gx1[83]*Gx2[183];
Gx3[72] = + Gx1[70]*Gx2[2] + Gx1[71]*Gx2[16] + Gx1[72]*Gx2[30] + Gx1[73]*Gx2[44] + Gx1[74]*Gx2[58] + Gx1[75]*Gx2[72] + Gx1[76]*Gx2[86] + Gx1[77]*Gx2[100] + Gx1[78]*Gx2[114] + Gx1[79]*Gx2[128] + Gx1[80]*Gx2[142] + Gx1[81]*Gx2[156] + Gx1[82]*Gx2[170] + Gx1[83]*Gx2[184];
Gx3[73] = + Gx1[70]*Gx2[3] + Gx1[71]*Gx2[17] + Gx1[72]*Gx2[31] + Gx1[73]*Gx2[45] + Gx1[74]*Gx2[59] + Gx1[75]*Gx2[73] + Gx1[76]*Gx2[87] + Gx1[77]*Gx2[101] + Gx1[78]*Gx2[115] + Gx1[79]*Gx2[129] + Gx1[80]*Gx2[143] + Gx1[81]*Gx2[157] + Gx1[82]*Gx2[171] + Gx1[83]*Gx2[185];
Gx3[74] = + Gx1[70]*Gx2[4] + Gx1[71]*Gx2[18] + Gx1[72]*Gx2[32] + Gx1[73]*Gx2[46] + Gx1[74]*Gx2[60] + Gx1[75]*Gx2[74] + Gx1[76]*Gx2[88] + Gx1[77]*Gx2[102] + Gx1[78]*Gx2[116] + Gx1[79]*Gx2[130] + Gx1[80]*Gx2[144] + Gx1[81]*Gx2[158] + Gx1[82]*Gx2[172] + Gx1[83]*Gx2[186];
Gx3[75] = + Gx1[70]*Gx2[5] + Gx1[71]*Gx2[19] + Gx1[72]*Gx2[33] + Gx1[73]*Gx2[47] + Gx1[74]*Gx2[61] + Gx1[75]*Gx2[75] + Gx1[76]*Gx2[89] + Gx1[77]*Gx2[103] + Gx1[78]*Gx2[117] + Gx1[79]*Gx2[131] + Gx1[80]*Gx2[145] + Gx1[81]*Gx2[159] + Gx1[82]*Gx2[173] + Gx1[83]*Gx2[187];
Gx3[76] = + Gx1[70]*Gx2[6] + Gx1[71]*Gx2[20] + Gx1[72]*Gx2[34] + Gx1[73]*Gx2[48] + Gx1[74]*Gx2[62] + Gx1[75]*Gx2[76] + Gx1[76]*Gx2[90] + Gx1[77]*Gx2[104] + Gx1[78]*Gx2[118] + Gx1[79]*Gx2[132] + Gx1[80]*Gx2[146] + Gx1[81]*Gx2[160] + Gx1[82]*Gx2[174] + Gx1[83]*Gx2[188];
Gx3[77] = + Gx1[70]*Gx2[7] + Gx1[71]*Gx2[21] + Gx1[72]*Gx2[35] + Gx1[73]*Gx2[49] + Gx1[74]*Gx2[63] + Gx1[75]*Gx2[77] + Gx1[76]*Gx2[91] + Gx1[77]*Gx2[105] + Gx1[78]*Gx2[119] + Gx1[79]*Gx2[133] + Gx1[80]*Gx2[147] + Gx1[81]*Gx2[161] + Gx1[82]*Gx2[175] + Gx1[83]*Gx2[189];
Gx3[78] = + Gx1[70]*Gx2[8] + Gx1[71]*Gx2[22] + Gx1[72]*Gx2[36] + Gx1[73]*Gx2[50] + Gx1[74]*Gx2[64] + Gx1[75]*Gx2[78] + Gx1[76]*Gx2[92] + Gx1[77]*Gx2[106] + Gx1[78]*Gx2[120] + Gx1[79]*Gx2[134] + Gx1[80]*Gx2[148] + Gx1[81]*Gx2[162] + Gx1[82]*Gx2[176] + Gx1[83]*Gx2[190];
Gx3[79] = + Gx1[70]*Gx2[9] + Gx1[71]*Gx2[23] + Gx1[72]*Gx2[37] + Gx1[73]*Gx2[51] + Gx1[74]*Gx2[65] + Gx1[75]*Gx2[79] + Gx1[76]*Gx2[93] + Gx1[77]*Gx2[107] + Gx1[78]*Gx2[121] + Gx1[79]*Gx2[135] + Gx1[80]*Gx2[149] + Gx1[81]*Gx2[163] + Gx1[82]*Gx2[177] + Gx1[83]*Gx2[191];
Gx3[80] = + Gx1[70]*Gx2[10] + Gx1[71]*Gx2[24] + Gx1[72]*Gx2[38] + Gx1[73]*Gx2[52] + Gx1[74]*Gx2[66] + Gx1[75]*Gx2[80] + Gx1[76]*Gx2[94] + Gx1[77]*Gx2[108] + Gx1[78]*Gx2[122] + Gx1[79]*Gx2[136] + Gx1[80]*Gx2[150] + Gx1[81]*Gx2[164] + Gx1[82]*Gx2[178] + Gx1[83]*Gx2[192];
Gx3[81] = + Gx1[70]*Gx2[11] + Gx1[71]*Gx2[25] + Gx1[72]*Gx2[39] + Gx1[73]*Gx2[53] + Gx1[74]*Gx2[67] + Gx1[75]*Gx2[81] + Gx1[76]*Gx2[95] + Gx1[77]*Gx2[109] + Gx1[78]*Gx2[123] + Gx1[79]*Gx2[137] + Gx1[80]*Gx2[151] + Gx1[81]*Gx2[165] + Gx1[82]*Gx2[179] + Gx1[83]*Gx2[193];
Gx3[82] = + Gx1[70]*Gx2[12] + Gx1[71]*Gx2[26] + Gx1[72]*Gx2[40] + Gx1[73]*Gx2[54] + Gx1[74]*Gx2[68] + Gx1[75]*Gx2[82] + Gx1[76]*Gx2[96] + Gx1[77]*Gx2[110] + Gx1[78]*Gx2[124] + Gx1[79]*Gx2[138] + Gx1[80]*Gx2[152] + Gx1[81]*Gx2[166] + Gx1[82]*Gx2[180] + Gx1[83]*Gx2[194];
Gx3[83] = + Gx1[70]*Gx2[13] + Gx1[71]*Gx2[27] + Gx1[72]*Gx2[41] + Gx1[73]*Gx2[55] + Gx1[74]*Gx2[69] + Gx1[75]*Gx2[83] + Gx1[76]*Gx2[97] + Gx1[77]*Gx2[111] + Gx1[78]*Gx2[125] + Gx1[79]*Gx2[139] + Gx1[80]*Gx2[153] + Gx1[81]*Gx2[167] + Gx1[82]*Gx2[181] + Gx1[83]*Gx2[195];
Gx3[84] = + Gx1[84]*Gx2[0] + Gx1[85]*Gx2[14] + Gx1[86]*Gx2[28] + Gx1[87]*Gx2[42] + Gx1[88]*Gx2[56] + Gx1[89]*Gx2[70] + Gx1[90]*Gx2[84] + Gx1[91]*Gx2[98] + Gx1[92]*Gx2[112] + Gx1[93]*Gx2[126] + Gx1[94]*Gx2[140] + Gx1[95]*Gx2[154] + Gx1[96]*Gx2[168] + Gx1[97]*Gx2[182];
Gx3[85] = + Gx1[84]*Gx2[1] + Gx1[85]*Gx2[15] + Gx1[86]*Gx2[29] + Gx1[87]*Gx2[43] + Gx1[88]*Gx2[57] + Gx1[89]*Gx2[71] + Gx1[90]*Gx2[85] + Gx1[91]*Gx2[99] + Gx1[92]*Gx2[113] + Gx1[93]*Gx2[127] + Gx1[94]*Gx2[141] + Gx1[95]*Gx2[155] + Gx1[96]*Gx2[169] + Gx1[97]*Gx2[183];
Gx3[86] = + Gx1[84]*Gx2[2] + Gx1[85]*Gx2[16] + Gx1[86]*Gx2[30] + Gx1[87]*Gx2[44] + Gx1[88]*Gx2[58] + Gx1[89]*Gx2[72] + Gx1[90]*Gx2[86] + Gx1[91]*Gx2[100] + Gx1[92]*Gx2[114] + Gx1[93]*Gx2[128] + Gx1[94]*Gx2[142] + Gx1[95]*Gx2[156] + Gx1[96]*Gx2[170] + Gx1[97]*Gx2[184];
Gx3[87] = + Gx1[84]*Gx2[3] + Gx1[85]*Gx2[17] + Gx1[86]*Gx2[31] + Gx1[87]*Gx2[45] + Gx1[88]*Gx2[59] + Gx1[89]*Gx2[73] + Gx1[90]*Gx2[87] + Gx1[91]*Gx2[101] + Gx1[92]*Gx2[115] + Gx1[93]*Gx2[129] + Gx1[94]*Gx2[143] + Gx1[95]*Gx2[157] + Gx1[96]*Gx2[171] + Gx1[97]*Gx2[185];
Gx3[88] = + Gx1[84]*Gx2[4] + Gx1[85]*Gx2[18] + Gx1[86]*Gx2[32] + Gx1[87]*Gx2[46] + Gx1[88]*Gx2[60] + Gx1[89]*Gx2[74] + Gx1[90]*Gx2[88] + Gx1[91]*Gx2[102] + Gx1[92]*Gx2[116] + Gx1[93]*Gx2[130] + Gx1[94]*Gx2[144] + Gx1[95]*Gx2[158] + Gx1[96]*Gx2[172] + Gx1[97]*Gx2[186];
Gx3[89] = + Gx1[84]*Gx2[5] + Gx1[85]*Gx2[19] + Gx1[86]*Gx2[33] + Gx1[87]*Gx2[47] + Gx1[88]*Gx2[61] + Gx1[89]*Gx2[75] + Gx1[90]*Gx2[89] + Gx1[91]*Gx2[103] + Gx1[92]*Gx2[117] + Gx1[93]*Gx2[131] + Gx1[94]*Gx2[145] + Gx1[95]*Gx2[159] + Gx1[96]*Gx2[173] + Gx1[97]*Gx2[187];
Gx3[90] = + Gx1[84]*Gx2[6] + Gx1[85]*Gx2[20] + Gx1[86]*Gx2[34] + Gx1[87]*Gx2[48] + Gx1[88]*Gx2[62] + Gx1[89]*Gx2[76] + Gx1[90]*Gx2[90] + Gx1[91]*Gx2[104] + Gx1[92]*Gx2[118] + Gx1[93]*Gx2[132] + Gx1[94]*Gx2[146] + Gx1[95]*Gx2[160] + Gx1[96]*Gx2[174] + Gx1[97]*Gx2[188];
Gx3[91] = + Gx1[84]*Gx2[7] + Gx1[85]*Gx2[21] + Gx1[86]*Gx2[35] + Gx1[87]*Gx2[49] + Gx1[88]*Gx2[63] + Gx1[89]*Gx2[77] + Gx1[90]*Gx2[91] + Gx1[91]*Gx2[105] + Gx1[92]*Gx2[119] + Gx1[93]*Gx2[133] + Gx1[94]*Gx2[147] + Gx1[95]*Gx2[161] + Gx1[96]*Gx2[175] + Gx1[97]*Gx2[189];
Gx3[92] = + Gx1[84]*Gx2[8] + Gx1[85]*Gx2[22] + Gx1[86]*Gx2[36] + Gx1[87]*Gx2[50] + Gx1[88]*Gx2[64] + Gx1[89]*Gx2[78] + Gx1[90]*Gx2[92] + Gx1[91]*Gx2[106] + Gx1[92]*Gx2[120] + Gx1[93]*Gx2[134] + Gx1[94]*Gx2[148] + Gx1[95]*Gx2[162] + Gx1[96]*Gx2[176] + Gx1[97]*Gx2[190];
Gx3[93] = + Gx1[84]*Gx2[9] + Gx1[85]*Gx2[23] + Gx1[86]*Gx2[37] + Gx1[87]*Gx2[51] + Gx1[88]*Gx2[65] + Gx1[89]*Gx2[79] + Gx1[90]*Gx2[93] + Gx1[91]*Gx2[107] + Gx1[92]*Gx2[121] + Gx1[93]*Gx2[135] + Gx1[94]*Gx2[149] + Gx1[95]*Gx2[163] + Gx1[96]*Gx2[177] + Gx1[97]*Gx2[191];
Gx3[94] = + Gx1[84]*Gx2[10] + Gx1[85]*Gx2[24] + Gx1[86]*Gx2[38] + Gx1[87]*Gx2[52] + Gx1[88]*Gx2[66] + Gx1[89]*Gx2[80] + Gx1[90]*Gx2[94] + Gx1[91]*Gx2[108] + Gx1[92]*Gx2[122] + Gx1[93]*Gx2[136] + Gx1[94]*Gx2[150] + Gx1[95]*Gx2[164] + Gx1[96]*Gx2[178] + Gx1[97]*Gx2[192];
Gx3[95] = + Gx1[84]*Gx2[11] + Gx1[85]*Gx2[25] + Gx1[86]*Gx2[39] + Gx1[87]*Gx2[53] + Gx1[88]*Gx2[67] + Gx1[89]*Gx2[81] + Gx1[90]*Gx2[95] + Gx1[91]*Gx2[109] + Gx1[92]*Gx2[123] + Gx1[93]*Gx2[137] + Gx1[94]*Gx2[151] + Gx1[95]*Gx2[165] + Gx1[96]*Gx2[179] + Gx1[97]*Gx2[193];
Gx3[96] = + Gx1[84]*Gx2[12] + Gx1[85]*Gx2[26] + Gx1[86]*Gx2[40] + Gx1[87]*Gx2[54] + Gx1[88]*Gx2[68] + Gx1[89]*Gx2[82] + Gx1[90]*Gx2[96] + Gx1[91]*Gx2[110] + Gx1[92]*Gx2[124] + Gx1[93]*Gx2[138] + Gx1[94]*Gx2[152] + Gx1[95]*Gx2[166] + Gx1[96]*Gx2[180] + Gx1[97]*Gx2[194];
Gx3[97] = + Gx1[84]*Gx2[13] + Gx1[85]*Gx2[27] + Gx1[86]*Gx2[41] + Gx1[87]*Gx2[55] + Gx1[88]*Gx2[69] + Gx1[89]*Gx2[83] + Gx1[90]*Gx2[97] + Gx1[91]*Gx2[111] + Gx1[92]*Gx2[125] + Gx1[93]*Gx2[139] + Gx1[94]*Gx2[153] + Gx1[95]*Gx2[167] + Gx1[96]*Gx2[181] + Gx1[97]*Gx2[195];
Gx3[98] = + Gx1[98]*Gx2[0] + Gx1[99]*Gx2[14] + Gx1[100]*Gx2[28] + Gx1[101]*Gx2[42] + Gx1[102]*Gx2[56] + Gx1[103]*Gx2[70] + Gx1[104]*Gx2[84] + Gx1[105]*Gx2[98] + Gx1[106]*Gx2[112] + Gx1[107]*Gx2[126] + Gx1[108]*Gx2[140] + Gx1[109]*Gx2[154] + Gx1[110]*Gx2[168] + Gx1[111]*Gx2[182];
Gx3[99] = + Gx1[98]*Gx2[1] + Gx1[99]*Gx2[15] + Gx1[100]*Gx2[29] + Gx1[101]*Gx2[43] + Gx1[102]*Gx2[57] + Gx1[103]*Gx2[71] + Gx1[104]*Gx2[85] + Gx1[105]*Gx2[99] + Gx1[106]*Gx2[113] + Gx1[107]*Gx2[127] + Gx1[108]*Gx2[141] + Gx1[109]*Gx2[155] + Gx1[110]*Gx2[169] + Gx1[111]*Gx2[183];
Gx3[100] = + Gx1[98]*Gx2[2] + Gx1[99]*Gx2[16] + Gx1[100]*Gx2[30] + Gx1[101]*Gx2[44] + Gx1[102]*Gx2[58] + Gx1[103]*Gx2[72] + Gx1[104]*Gx2[86] + Gx1[105]*Gx2[100] + Gx1[106]*Gx2[114] + Gx1[107]*Gx2[128] + Gx1[108]*Gx2[142] + Gx1[109]*Gx2[156] + Gx1[110]*Gx2[170] + Gx1[111]*Gx2[184];
Gx3[101] = + Gx1[98]*Gx2[3] + Gx1[99]*Gx2[17] + Gx1[100]*Gx2[31] + Gx1[101]*Gx2[45] + Gx1[102]*Gx2[59] + Gx1[103]*Gx2[73] + Gx1[104]*Gx2[87] + Gx1[105]*Gx2[101] + Gx1[106]*Gx2[115] + Gx1[107]*Gx2[129] + Gx1[108]*Gx2[143] + Gx1[109]*Gx2[157] + Gx1[110]*Gx2[171] + Gx1[111]*Gx2[185];
Gx3[102] = + Gx1[98]*Gx2[4] + Gx1[99]*Gx2[18] + Gx1[100]*Gx2[32] + Gx1[101]*Gx2[46] + Gx1[102]*Gx2[60] + Gx1[103]*Gx2[74] + Gx1[104]*Gx2[88] + Gx1[105]*Gx2[102] + Gx1[106]*Gx2[116] + Gx1[107]*Gx2[130] + Gx1[108]*Gx2[144] + Gx1[109]*Gx2[158] + Gx1[110]*Gx2[172] + Gx1[111]*Gx2[186];
Gx3[103] = + Gx1[98]*Gx2[5] + Gx1[99]*Gx2[19] + Gx1[100]*Gx2[33] + Gx1[101]*Gx2[47] + Gx1[102]*Gx2[61] + Gx1[103]*Gx2[75] + Gx1[104]*Gx2[89] + Gx1[105]*Gx2[103] + Gx1[106]*Gx2[117] + Gx1[107]*Gx2[131] + Gx1[108]*Gx2[145] + Gx1[109]*Gx2[159] + Gx1[110]*Gx2[173] + Gx1[111]*Gx2[187];
Gx3[104] = + Gx1[98]*Gx2[6] + Gx1[99]*Gx2[20] + Gx1[100]*Gx2[34] + Gx1[101]*Gx2[48] + Gx1[102]*Gx2[62] + Gx1[103]*Gx2[76] + Gx1[104]*Gx2[90] + Gx1[105]*Gx2[104] + Gx1[106]*Gx2[118] + Gx1[107]*Gx2[132] + Gx1[108]*Gx2[146] + Gx1[109]*Gx2[160] + Gx1[110]*Gx2[174] + Gx1[111]*Gx2[188];
Gx3[105] = + Gx1[98]*Gx2[7] + Gx1[99]*Gx2[21] + Gx1[100]*Gx2[35] + Gx1[101]*Gx2[49] + Gx1[102]*Gx2[63] + Gx1[103]*Gx2[77] + Gx1[104]*Gx2[91] + Gx1[105]*Gx2[105] + Gx1[106]*Gx2[119] + Gx1[107]*Gx2[133] + Gx1[108]*Gx2[147] + Gx1[109]*Gx2[161] + Gx1[110]*Gx2[175] + Gx1[111]*Gx2[189];
Gx3[106] = + Gx1[98]*Gx2[8] + Gx1[99]*Gx2[22] + Gx1[100]*Gx2[36] + Gx1[101]*Gx2[50] + Gx1[102]*Gx2[64] + Gx1[103]*Gx2[78] + Gx1[104]*Gx2[92] + Gx1[105]*Gx2[106] + Gx1[106]*Gx2[120] + Gx1[107]*Gx2[134] + Gx1[108]*Gx2[148] + Gx1[109]*Gx2[162] + Gx1[110]*Gx2[176] + Gx1[111]*Gx2[190];
Gx3[107] = + Gx1[98]*Gx2[9] + Gx1[99]*Gx2[23] + Gx1[100]*Gx2[37] + Gx1[101]*Gx2[51] + Gx1[102]*Gx2[65] + Gx1[103]*Gx2[79] + Gx1[104]*Gx2[93] + Gx1[105]*Gx2[107] + Gx1[106]*Gx2[121] + Gx1[107]*Gx2[135] + Gx1[108]*Gx2[149] + Gx1[109]*Gx2[163] + Gx1[110]*Gx2[177] + Gx1[111]*Gx2[191];
Gx3[108] = + Gx1[98]*Gx2[10] + Gx1[99]*Gx2[24] + Gx1[100]*Gx2[38] + Gx1[101]*Gx2[52] + Gx1[102]*Gx2[66] + Gx1[103]*Gx2[80] + Gx1[104]*Gx2[94] + Gx1[105]*Gx2[108] + Gx1[106]*Gx2[122] + Gx1[107]*Gx2[136] + Gx1[108]*Gx2[150] + Gx1[109]*Gx2[164] + Gx1[110]*Gx2[178] + Gx1[111]*Gx2[192];
Gx3[109] = + Gx1[98]*Gx2[11] + Gx1[99]*Gx2[25] + Gx1[100]*Gx2[39] + Gx1[101]*Gx2[53] + Gx1[102]*Gx2[67] + Gx1[103]*Gx2[81] + Gx1[104]*Gx2[95] + Gx1[105]*Gx2[109] + Gx1[106]*Gx2[123] + Gx1[107]*Gx2[137] + Gx1[108]*Gx2[151] + Gx1[109]*Gx2[165] + Gx1[110]*Gx2[179] + Gx1[111]*Gx2[193];
Gx3[110] = + Gx1[98]*Gx2[12] + Gx1[99]*Gx2[26] + Gx1[100]*Gx2[40] + Gx1[101]*Gx2[54] + Gx1[102]*Gx2[68] + Gx1[103]*Gx2[82] + Gx1[104]*Gx2[96] + Gx1[105]*Gx2[110] + Gx1[106]*Gx2[124] + Gx1[107]*Gx2[138] + Gx1[108]*Gx2[152] + Gx1[109]*Gx2[166] + Gx1[110]*Gx2[180] + Gx1[111]*Gx2[194];
Gx3[111] = + Gx1[98]*Gx2[13] + Gx1[99]*Gx2[27] + Gx1[100]*Gx2[41] + Gx1[101]*Gx2[55] + Gx1[102]*Gx2[69] + Gx1[103]*Gx2[83] + Gx1[104]*Gx2[97] + Gx1[105]*Gx2[111] + Gx1[106]*Gx2[125] + Gx1[107]*Gx2[139] + Gx1[108]*Gx2[153] + Gx1[109]*Gx2[167] + Gx1[110]*Gx2[181] + Gx1[111]*Gx2[195];
Gx3[112] = + Gx1[112]*Gx2[0] + Gx1[113]*Gx2[14] + Gx1[114]*Gx2[28] + Gx1[115]*Gx2[42] + Gx1[116]*Gx2[56] + Gx1[117]*Gx2[70] + Gx1[118]*Gx2[84] + Gx1[119]*Gx2[98] + Gx1[120]*Gx2[112] + Gx1[121]*Gx2[126] + Gx1[122]*Gx2[140] + Gx1[123]*Gx2[154] + Gx1[124]*Gx2[168] + Gx1[125]*Gx2[182];
Gx3[113] = + Gx1[112]*Gx2[1] + Gx1[113]*Gx2[15] + Gx1[114]*Gx2[29] + Gx1[115]*Gx2[43] + Gx1[116]*Gx2[57] + Gx1[117]*Gx2[71] + Gx1[118]*Gx2[85] + Gx1[119]*Gx2[99] + Gx1[120]*Gx2[113] + Gx1[121]*Gx2[127] + Gx1[122]*Gx2[141] + Gx1[123]*Gx2[155] + Gx1[124]*Gx2[169] + Gx1[125]*Gx2[183];
Gx3[114] = + Gx1[112]*Gx2[2] + Gx1[113]*Gx2[16] + Gx1[114]*Gx2[30] + Gx1[115]*Gx2[44] + Gx1[116]*Gx2[58] + Gx1[117]*Gx2[72] + Gx1[118]*Gx2[86] + Gx1[119]*Gx2[100] + Gx1[120]*Gx2[114] + Gx1[121]*Gx2[128] + Gx1[122]*Gx2[142] + Gx1[123]*Gx2[156] + Gx1[124]*Gx2[170] + Gx1[125]*Gx2[184];
Gx3[115] = + Gx1[112]*Gx2[3] + Gx1[113]*Gx2[17] + Gx1[114]*Gx2[31] + Gx1[115]*Gx2[45] + Gx1[116]*Gx2[59] + Gx1[117]*Gx2[73] + Gx1[118]*Gx2[87] + Gx1[119]*Gx2[101] + Gx1[120]*Gx2[115] + Gx1[121]*Gx2[129] + Gx1[122]*Gx2[143] + Gx1[123]*Gx2[157] + Gx1[124]*Gx2[171] + Gx1[125]*Gx2[185];
Gx3[116] = + Gx1[112]*Gx2[4] + Gx1[113]*Gx2[18] + Gx1[114]*Gx2[32] + Gx1[115]*Gx2[46] + Gx1[116]*Gx2[60] + Gx1[117]*Gx2[74] + Gx1[118]*Gx2[88] + Gx1[119]*Gx2[102] + Gx1[120]*Gx2[116] + Gx1[121]*Gx2[130] + Gx1[122]*Gx2[144] + Gx1[123]*Gx2[158] + Gx1[124]*Gx2[172] + Gx1[125]*Gx2[186];
Gx3[117] = + Gx1[112]*Gx2[5] + Gx1[113]*Gx2[19] + Gx1[114]*Gx2[33] + Gx1[115]*Gx2[47] + Gx1[116]*Gx2[61] + Gx1[117]*Gx2[75] + Gx1[118]*Gx2[89] + Gx1[119]*Gx2[103] + Gx1[120]*Gx2[117] + Gx1[121]*Gx2[131] + Gx1[122]*Gx2[145] + Gx1[123]*Gx2[159] + Gx1[124]*Gx2[173] + Gx1[125]*Gx2[187];
Gx3[118] = + Gx1[112]*Gx2[6] + Gx1[113]*Gx2[20] + Gx1[114]*Gx2[34] + Gx1[115]*Gx2[48] + Gx1[116]*Gx2[62] + Gx1[117]*Gx2[76] + Gx1[118]*Gx2[90] + Gx1[119]*Gx2[104] + Gx1[120]*Gx2[118] + Gx1[121]*Gx2[132] + Gx1[122]*Gx2[146] + Gx1[123]*Gx2[160] + Gx1[124]*Gx2[174] + Gx1[125]*Gx2[188];
Gx3[119] = + Gx1[112]*Gx2[7] + Gx1[113]*Gx2[21] + Gx1[114]*Gx2[35] + Gx1[115]*Gx2[49] + Gx1[116]*Gx2[63] + Gx1[117]*Gx2[77] + Gx1[118]*Gx2[91] + Gx1[119]*Gx2[105] + Gx1[120]*Gx2[119] + Gx1[121]*Gx2[133] + Gx1[122]*Gx2[147] + Gx1[123]*Gx2[161] + Gx1[124]*Gx2[175] + Gx1[125]*Gx2[189];
Gx3[120] = + Gx1[112]*Gx2[8] + Gx1[113]*Gx2[22] + Gx1[114]*Gx2[36] + Gx1[115]*Gx2[50] + Gx1[116]*Gx2[64] + Gx1[117]*Gx2[78] + Gx1[118]*Gx2[92] + Gx1[119]*Gx2[106] + Gx1[120]*Gx2[120] + Gx1[121]*Gx2[134] + Gx1[122]*Gx2[148] + Gx1[123]*Gx2[162] + Gx1[124]*Gx2[176] + Gx1[125]*Gx2[190];
Gx3[121] = + Gx1[112]*Gx2[9] + Gx1[113]*Gx2[23] + Gx1[114]*Gx2[37] + Gx1[115]*Gx2[51] + Gx1[116]*Gx2[65] + Gx1[117]*Gx2[79] + Gx1[118]*Gx2[93] + Gx1[119]*Gx2[107] + Gx1[120]*Gx2[121] + Gx1[121]*Gx2[135] + Gx1[122]*Gx2[149] + Gx1[123]*Gx2[163] + Gx1[124]*Gx2[177] + Gx1[125]*Gx2[191];
Gx3[122] = + Gx1[112]*Gx2[10] + Gx1[113]*Gx2[24] + Gx1[114]*Gx2[38] + Gx1[115]*Gx2[52] + Gx1[116]*Gx2[66] + Gx1[117]*Gx2[80] + Gx1[118]*Gx2[94] + Gx1[119]*Gx2[108] + Gx1[120]*Gx2[122] + Gx1[121]*Gx2[136] + Gx1[122]*Gx2[150] + Gx1[123]*Gx2[164] + Gx1[124]*Gx2[178] + Gx1[125]*Gx2[192];
Gx3[123] = + Gx1[112]*Gx2[11] + Gx1[113]*Gx2[25] + Gx1[114]*Gx2[39] + Gx1[115]*Gx2[53] + Gx1[116]*Gx2[67] + Gx1[117]*Gx2[81] + Gx1[118]*Gx2[95] + Gx1[119]*Gx2[109] + Gx1[120]*Gx2[123] + Gx1[121]*Gx2[137] + Gx1[122]*Gx2[151] + Gx1[123]*Gx2[165] + Gx1[124]*Gx2[179] + Gx1[125]*Gx2[193];
Gx3[124] = + Gx1[112]*Gx2[12] + Gx1[113]*Gx2[26] + Gx1[114]*Gx2[40] + Gx1[115]*Gx2[54] + Gx1[116]*Gx2[68] + Gx1[117]*Gx2[82] + Gx1[118]*Gx2[96] + Gx1[119]*Gx2[110] + Gx1[120]*Gx2[124] + Gx1[121]*Gx2[138] + Gx1[122]*Gx2[152] + Gx1[123]*Gx2[166] + Gx1[124]*Gx2[180] + Gx1[125]*Gx2[194];
Gx3[125] = + Gx1[112]*Gx2[13] + Gx1[113]*Gx2[27] + Gx1[114]*Gx2[41] + Gx1[115]*Gx2[55] + Gx1[116]*Gx2[69] + Gx1[117]*Gx2[83] + Gx1[118]*Gx2[97] + Gx1[119]*Gx2[111] + Gx1[120]*Gx2[125] + Gx1[121]*Gx2[139] + Gx1[122]*Gx2[153] + Gx1[123]*Gx2[167] + Gx1[124]*Gx2[181] + Gx1[125]*Gx2[195];
Gx3[126] = + Gx1[126]*Gx2[0] + Gx1[127]*Gx2[14] + Gx1[128]*Gx2[28] + Gx1[129]*Gx2[42] + Gx1[130]*Gx2[56] + Gx1[131]*Gx2[70] + Gx1[132]*Gx2[84] + Gx1[133]*Gx2[98] + Gx1[134]*Gx2[112] + Gx1[135]*Gx2[126] + Gx1[136]*Gx2[140] + Gx1[137]*Gx2[154] + Gx1[138]*Gx2[168] + Gx1[139]*Gx2[182];
Gx3[127] = + Gx1[126]*Gx2[1] + Gx1[127]*Gx2[15] + Gx1[128]*Gx2[29] + Gx1[129]*Gx2[43] + Gx1[130]*Gx2[57] + Gx1[131]*Gx2[71] + Gx1[132]*Gx2[85] + Gx1[133]*Gx2[99] + Gx1[134]*Gx2[113] + Gx1[135]*Gx2[127] + Gx1[136]*Gx2[141] + Gx1[137]*Gx2[155] + Gx1[138]*Gx2[169] + Gx1[139]*Gx2[183];
Gx3[128] = + Gx1[126]*Gx2[2] + Gx1[127]*Gx2[16] + Gx1[128]*Gx2[30] + Gx1[129]*Gx2[44] + Gx1[130]*Gx2[58] + Gx1[131]*Gx2[72] + Gx1[132]*Gx2[86] + Gx1[133]*Gx2[100] + Gx1[134]*Gx2[114] + Gx1[135]*Gx2[128] + Gx1[136]*Gx2[142] + Gx1[137]*Gx2[156] + Gx1[138]*Gx2[170] + Gx1[139]*Gx2[184];
Gx3[129] = + Gx1[126]*Gx2[3] + Gx1[127]*Gx2[17] + Gx1[128]*Gx2[31] + Gx1[129]*Gx2[45] + Gx1[130]*Gx2[59] + Gx1[131]*Gx2[73] + Gx1[132]*Gx2[87] + Gx1[133]*Gx2[101] + Gx1[134]*Gx2[115] + Gx1[135]*Gx2[129] + Gx1[136]*Gx2[143] + Gx1[137]*Gx2[157] + Gx1[138]*Gx2[171] + Gx1[139]*Gx2[185];
Gx3[130] = + Gx1[126]*Gx2[4] + Gx1[127]*Gx2[18] + Gx1[128]*Gx2[32] + Gx1[129]*Gx2[46] + Gx1[130]*Gx2[60] + Gx1[131]*Gx2[74] + Gx1[132]*Gx2[88] + Gx1[133]*Gx2[102] + Gx1[134]*Gx2[116] + Gx1[135]*Gx2[130] + Gx1[136]*Gx2[144] + Gx1[137]*Gx2[158] + Gx1[138]*Gx2[172] + Gx1[139]*Gx2[186];
Gx3[131] = + Gx1[126]*Gx2[5] + Gx1[127]*Gx2[19] + Gx1[128]*Gx2[33] + Gx1[129]*Gx2[47] + Gx1[130]*Gx2[61] + Gx1[131]*Gx2[75] + Gx1[132]*Gx2[89] + Gx1[133]*Gx2[103] + Gx1[134]*Gx2[117] + Gx1[135]*Gx2[131] + Gx1[136]*Gx2[145] + Gx1[137]*Gx2[159] + Gx1[138]*Gx2[173] + Gx1[139]*Gx2[187];
Gx3[132] = + Gx1[126]*Gx2[6] + Gx1[127]*Gx2[20] + Gx1[128]*Gx2[34] + Gx1[129]*Gx2[48] + Gx1[130]*Gx2[62] + Gx1[131]*Gx2[76] + Gx1[132]*Gx2[90] + Gx1[133]*Gx2[104] + Gx1[134]*Gx2[118] + Gx1[135]*Gx2[132] + Gx1[136]*Gx2[146] + Gx1[137]*Gx2[160] + Gx1[138]*Gx2[174] + Gx1[139]*Gx2[188];
Gx3[133] = + Gx1[126]*Gx2[7] + Gx1[127]*Gx2[21] + Gx1[128]*Gx2[35] + Gx1[129]*Gx2[49] + Gx1[130]*Gx2[63] + Gx1[131]*Gx2[77] + Gx1[132]*Gx2[91] + Gx1[133]*Gx2[105] + Gx1[134]*Gx2[119] + Gx1[135]*Gx2[133] + Gx1[136]*Gx2[147] + Gx1[137]*Gx2[161] + Gx1[138]*Gx2[175] + Gx1[139]*Gx2[189];
Gx3[134] = + Gx1[126]*Gx2[8] + Gx1[127]*Gx2[22] + Gx1[128]*Gx2[36] + Gx1[129]*Gx2[50] + Gx1[130]*Gx2[64] + Gx1[131]*Gx2[78] + Gx1[132]*Gx2[92] + Gx1[133]*Gx2[106] + Gx1[134]*Gx2[120] + Gx1[135]*Gx2[134] + Gx1[136]*Gx2[148] + Gx1[137]*Gx2[162] + Gx1[138]*Gx2[176] + Gx1[139]*Gx2[190];
Gx3[135] = + Gx1[126]*Gx2[9] + Gx1[127]*Gx2[23] + Gx1[128]*Gx2[37] + Gx1[129]*Gx2[51] + Gx1[130]*Gx2[65] + Gx1[131]*Gx2[79] + Gx1[132]*Gx2[93] + Gx1[133]*Gx2[107] + Gx1[134]*Gx2[121] + Gx1[135]*Gx2[135] + Gx1[136]*Gx2[149] + Gx1[137]*Gx2[163] + Gx1[138]*Gx2[177] + Gx1[139]*Gx2[191];
Gx3[136] = + Gx1[126]*Gx2[10] + Gx1[127]*Gx2[24] + Gx1[128]*Gx2[38] + Gx1[129]*Gx2[52] + Gx1[130]*Gx2[66] + Gx1[131]*Gx2[80] + Gx1[132]*Gx2[94] + Gx1[133]*Gx2[108] + Gx1[134]*Gx2[122] + Gx1[135]*Gx2[136] + Gx1[136]*Gx2[150] + Gx1[137]*Gx2[164] + Gx1[138]*Gx2[178] + Gx1[139]*Gx2[192];
Gx3[137] = + Gx1[126]*Gx2[11] + Gx1[127]*Gx2[25] + Gx1[128]*Gx2[39] + Gx1[129]*Gx2[53] + Gx1[130]*Gx2[67] + Gx1[131]*Gx2[81] + Gx1[132]*Gx2[95] + Gx1[133]*Gx2[109] + Gx1[134]*Gx2[123] + Gx1[135]*Gx2[137] + Gx1[136]*Gx2[151] + Gx1[137]*Gx2[165] + Gx1[138]*Gx2[179] + Gx1[139]*Gx2[193];
Gx3[138] = + Gx1[126]*Gx2[12] + Gx1[127]*Gx2[26] + Gx1[128]*Gx2[40] + Gx1[129]*Gx2[54] + Gx1[130]*Gx2[68] + Gx1[131]*Gx2[82] + Gx1[132]*Gx2[96] + Gx1[133]*Gx2[110] + Gx1[134]*Gx2[124] + Gx1[135]*Gx2[138] + Gx1[136]*Gx2[152] + Gx1[137]*Gx2[166] + Gx1[138]*Gx2[180] + Gx1[139]*Gx2[194];
Gx3[139] = + Gx1[126]*Gx2[13] + Gx1[127]*Gx2[27] + Gx1[128]*Gx2[41] + Gx1[129]*Gx2[55] + Gx1[130]*Gx2[69] + Gx1[131]*Gx2[83] + Gx1[132]*Gx2[97] + Gx1[133]*Gx2[111] + Gx1[134]*Gx2[125] + Gx1[135]*Gx2[139] + Gx1[136]*Gx2[153] + Gx1[137]*Gx2[167] + Gx1[138]*Gx2[181] + Gx1[139]*Gx2[195];
Gx3[140] = + Gx1[140]*Gx2[0] + Gx1[141]*Gx2[14] + Gx1[142]*Gx2[28] + Gx1[143]*Gx2[42] + Gx1[144]*Gx2[56] + Gx1[145]*Gx2[70] + Gx1[146]*Gx2[84] + Gx1[147]*Gx2[98] + Gx1[148]*Gx2[112] + Gx1[149]*Gx2[126] + Gx1[150]*Gx2[140] + Gx1[151]*Gx2[154] + Gx1[152]*Gx2[168] + Gx1[153]*Gx2[182];
Gx3[141] = + Gx1[140]*Gx2[1] + Gx1[141]*Gx2[15] + Gx1[142]*Gx2[29] + Gx1[143]*Gx2[43] + Gx1[144]*Gx2[57] + Gx1[145]*Gx2[71] + Gx1[146]*Gx2[85] + Gx1[147]*Gx2[99] + Gx1[148]*Gx2[113] + Gx1[149]*Gx2[127] + Gx1[150]*Gx2[141] + Gx1[151]*Gx2[155] + Gx1[152]*Gx2[169] + Gx1[153]*Gx2[183];
Gx3[142] = + Gx1[140]*Gx2[2] + Gx1[141]*Gx2[16] + Gx1[142]*Gx2[30] + Gx1[143]*Gx2[44] + Gx1[144]*Gx2[58] + Gx1[145]*Gx2[72] + Gx1[146]*Gx2[86] + Gx1[147]*Gx2[100] + Gx1[148]*Gx2[114] + Gx1[149]*Gx2[128] + Gx1[150]*Gx2[142] + Gx1[151]*Gx2[156] + Gx1[152]*Gx2[170] + Gx1[153]*Gx2[184];
Gx3[143] = + Gx1[140]*Gx2[3] + Gx1[141]*Gx2[17] + Gx1[142]*Gx2[31] + Gx1[143]*Gx2[45] + Gx1[144]*Gx2[59] + Gx1[145]*Gx2[73] + Gx1[146]*Gx2[87] + Gx1[147]*Gx2[101] + Gx1[148]*Gx2[115] + Gx1[149]*Gx2[129] + Gx1[150]*Gx2[143] + Gx1[151]*Gx2[157] + Gx1[152]*Gx2[171] + Gx1[153]*Gx2[185];
Gx3[144] = + Gx1[140]*Gx2[4] + Gx1[141]*Gx2[18] + Gx1[142]*Gx2[32] + Gx1[143]*Gx2[46] + Gx1[144]*Gx2[60] + Gx1[145]*Gx2[74] + Gx1[146]*Gx2[88] + Gx1[147]*Gx2[102] + Gx1[148]*Gx2[116] + Gx1[149]*Gx2[130] + Gx1[150]*Gx2[144] + Gx1[151]*Gx2[158] + Gx1[152]*Gx2[172] + Gx1[153]*Gx2[186];
Gx3[145] = + Gx1[140]*Gx2[5] + Gx1[141]*Gx2[19] + Gx1[142]*Gx2[33] + Gx1[143]*Gx2[47] + Gx1[144]*Gx2[61] + Gx1[145]*Gx2[75] + Gx1[146]*Gx2[89] + Gx1[147]*Gx2[103] + Gx1[148]*Gx2[117] + Gx1[149]*Gx2[131] + Gx1[150]*Gx2[145] + Gx1[151]*Gx2[159] + Gx1[152]*Gx2[173] + Gx1[153]*Gx2[187];
Gx3[146] = + Gx1[140]*Gx2[6] + Gx1[141]*Gx2[20] + Gx1[142]*Gx2[34] + Gx1[143]*Gx2[48] + Gx1[144]*Gx2[62] + Gx1[145]*Gx2[76] + Gx1[146]*Gx2[90] + Gx1[147]*Gx2[104] + Gx1[148]*Gx2[118] + Gx1[149]*Gx2[132] + Gx1[150]*Gx2[146] + Gx1[151]*Gx2[160] + Gx1[152]*Gx2[174] + Gx1[153]*Gx2[188];
Gx3[147] = + Gx1[140]*Gx2[7] + Gx1[141]*Gx2[21] + Gx1[142]*Gx2[35] + Gx1[143]*Gx2[49] + Gx1[144]*Gx2[63] + Gx1[145]*Gx2[77] + Gx1[146]*Gx2[91] + Gx1[147]*Gx2[105] + Gx1[148]*Gx2[119] + Gx1[149]*Gx2[133] + Gx1[150]*Gx2[147] + Gx1[151]*Gx2[161] + Gx1[152]*Gx2[175] + Gx1[153]*Gx2[189];
Gx3[148] = + Gx1[140]*Gx2[8] + Gx1[141]*Gx2[22] + Gx1[142]*Gx2[36] + Gx1[143]*Gx2[50] + Gx1[144]*Gx2[64] + Gx1[145]*Gx2[78] + Gx1[146]*Gx2[92] + Gx1[147]*Gx2[106] + Gx1[148]*Gx2[120] + Gx1[149]*Gx2[134] + Gx1[150]*Gx2[148] + Gx1[151]*Gx2[162] + Gx1[152]*Gx2[176] + Gx1[153]*Gx2[190];
Gx3[149] = + Gx1[140]*Gx2[9] + Gx1[141]*Gx2[23] + Gx1[142]*Gx2[37] + Gx1[143]*Gx2[51] + Gx1[144]*Gx2[65] + Gx1[145]*Gx2[79] + Gx1[146]*Gx2[93] + Gx1[147]*Gx2[107] + Gx1[148]*Gx2[121] + Gx1[149]*Gx2[135] + Gx1[150]*Gx2[149] + Gx1[151]*Gx2[163] + Gx1[152]*Gx2[177] + Gx1[153]*Gx2[191];
Gx3[150] = + Gx1[140]*Gx2[10] + Gx1[141]*Gx2[24] + Gx1[142]*Gx2[38] + Gx1[143]*Gx2[52] + Gx1[144]*Gx2[66] + Gx1[145]*Gx2[80] + Gx1[146]*Gx2[94] + Gx1[147]*Gx2[108] + Gx1[148]*Gx2[122] + Gx1[149]*Gx2[136] + Gx1[150]*Gx2[150] + Gx1[151]*Gx2[164] + Gx1[152]*Gx2[178] + Gx1[153]*Gx2[192];
Gx3[151] = + Gx1[140]*Gx2[11] + Gx1[141]*Gx2[25] + Gx1[142]*Gx2[39] + Gx1[143]*Gx2[53] + Gx1[144]*Gx2[67] + Gx1[145]*Gx2[81] + Gx1[146]*Gx2[95] + Gx1[147]*Gx2[109] + Gx1[148]*Gx2[123] + Gx1[149]*Gx2[137] + Gx1[150]*Gx2[151] + Gx1[151]*Gx2[165] + Gx1[152]*Gx2[179] + Gx1[153]*Gx2[193];
Gx3[152] = + Gx1[140]*Gx2[12] + Gx1[141]*Gx2[26] + Gx1[142]*Gx2[40] + Gx1[143]*Gx2[54] + Gx1[144]*Gx2[68] + Gx1[145]*Gx2[82] + Gx1[146]*Gx2[96] + Gx1[147]*Gx2[110] + Gx1[148]*Gx2[124] + Gx1[149]*Gx2[138] + Gx1[150]*Gx2[152] + Gx1[151]*Gx2[166] + Gx1[152]*Gx2[180] + Gx1[153]*Gx2[194];
Gx3[153] = + Gx1[140]*Gx2[13] + Gx1[141]*Gx2[27] + Gx1[142]*Gx2[41] + Gx1[143]*Gx2[55] + Gx1[144]*Gx2[69] + Gx1[145]*Gx2[83] + Gx1[146]*Gx2[97] + Gx1[147]*Gx2[111] + Gx1[148]*Gx2[125] + Gx1[149]*Gx2[139] + Gx1[150]*Gx2[153] + Gx1[151]*Gx2[167] + Gx1[152]*Gx2[181] + Gx1[153]*Gx2[195];
Gx3[154] = + Gx1[154]*Gx2[0] + Gx1[155]*Gx2[14] + Gx1[156]*Gx2[28] + Gx1[157]*Gx2[42] + Gx1[158]*Gx2[56] + Gx1[159]*Gx2[70] + Gx1[160]*Gx2[84] + Gx1[161]*Gx2[98] + Gx1[162]*Gx2[112] + Gx1[163]*Gx2[126] + Gx1[164]*Gx2[140] + Gx1[165]*Gx2[154] + Gx1[166]*Gx2[168] + Gx1[167]*Gx2[182];
Gx3[155] = + Gx1[154]*Gx2[1] + Gx1[155]*Gx2[15] + Gx1[156]*Gx2[29] + Gx1[157]*Gx2[43] + Gx1[158]*Gx2[57] + Gx1[159]*Gx2[71] + Gx1[160]*Gx2[85] + Gx1[161]*Gx2[99] + Gx1[162]*Gx2[113] + Gx1[163]*Gx2[127] + Gx1[164]*Gx2[141] + Gx1[165]*Gx2[155] + Gx1[166]*Gx2[169] + Gx1[167]*Gx2[183];
Gx3[156] = + Gx1[154]*Gx2[2] + Gx1[155]*Gx2[16] + Gx1[156]*Gx2[30] + Gx1[157]*Gx2[44] + Gx1[158]*Gx2[58] + Gx1[159]*Gx2[72] + Gx1[160]*Gx2[86] + Gx1[161]*Gx2[100] + Gx1[162]*Gx2[114] + Gx1[163]*Gx2[128] + Gx1[164]*Gx2[142] + Gx1[165]*Gx2[156] + Gx1[166]*Gx2[170] + Gx1[167]*Gx2[184];
Gx3[157] = + Gx1[154]*Gx2[3] + Gx1[155]*Gx2[17] + Gx1[156]*Gx2[31] + Gx1[157]*Gx2[45] + Gx1[158]*Gx2[59] + Gx1[159]*Gx2[73] + Gx1[160]*Gx2[87] + Gx1[161]*Gx2[101] + Gx1[162]*Gx2[115] + Gx1[163]*Gx2[129] + Gx1[164]*Gx2[143] + Gx1[165]*Gx2[157] + Gx1[166]*Gx2[171] + Gx1[167]*Gx2[185];
Gx3[158] = + Gx1[154]*Gx2[4] + Gx1[155]*Gx2[18] + Gx1[156]*Gx2[32] + Gx1[157]*Gx2[46] + Gx1[158]*Gx2[60] + Gx1[159]*Gx2[74] + Gx1[160]*Gx2[88] + Gx1[161]*Gx2[102] + Gx1[162]*Gx2[116] + Gx1[163]*Gx2[130] + Gx1[164]*Gx2[144] + Gx1[165]*Gx2[158] + Gx1[166]*Gx2[172] + Gx1[167]*Gx2[186];
Gx3[159] = + Gx1[154]*Gx2[5] + Gx1[155]*Gx2[19] + Gx1[156]*Gx2[33] + Gx1[157]*Gx2[47] + Gx1[158]*Gx2[61] + Gx1[159]*Gx2[75] + Gx1[160]*Gx2[89] + Gx1[161]*Gx2[103] + Gx1[162]*Gx2[117] + Gx1[163]*Gx2[131] + Gx1[164]*Gx2[145] + Gx1[165]*Gx2[159] + Gx1[166]*Gx2[173] + Gx1[167]*Gx2[187];
Gx3[160] = + Gx1[154]*Gx2[6] + Gx1[155]*Gx2[20] + Gx1[156]*Gx2[34] + Gx1[157]*Gx2[48] + Gx1[158]*Gx2[62] + Gx1[159]*Gx2[76] + Gx1[160]*Gx2[90] + Gx1[161]*Gx2[104] + Gx1[162]*Gx2[118] + Gx1[163]*Gx2[132] + Gx1[164]*Gx2[146] + Gx1[165]*Gx2[160] + Gx1[166]*Gx2[174] + Gx1[167]*Gx2[188];
Gx3[161] = + Gx1[154]*Gx2[7] + Gx1[155]*Gx2[21] + Gx1[156]*Gx2[35] + Gx1[157]*Gx2[49] + Gx1[158]*Gx2[63] + Gx1[159]*Gx2[77] + Gx1[160]*Gx2[91] + Gx1[161]*Gx2[105] + Gx1[162]*Gx2[119] + Gx1[163]*Gx2[133] + Gx1[164]*Gx2[147] + Gx1[165]*Gx2[161] + Gx1[166]*Gx2[175] + Gx1[167]*Gx2[189];
Gx3[162] = + Gx1[154]*Gx2[8] + Gx1[155]*Gx2[22] + Gx1[156]*Gx2[36] + Gx1[157]*Gx2[50] + Gx1[158]*Gx2[64] + Gx1[159]*Gx2[78] + Gx1[160]*Gx2[92] + Gx1[161]*Gx2[106] + Gx1[162]*Gx2[120] + Gx1[163]*Gx2[134] + Gx1[164]*Gx2[148] + Gx1[165]*Gx2[162] + Gx1[166]*Gx2[176] + Gx1[167]*Gx2[190];
Gx3[163] = + Gx1[154]*Gx2[9] + Gx1[155]*Gx2[23] + Gx1[156]*Gx2[37] + Gx1[157]*Gx2[51] + Gx1[158]*Gx2[65] + Gx1[159]*Gx2[79] + Gx1[160]*Gx2[93] + Gx1[161]*Gx2[107] + Gx1[162]*Gx2[121] + Gx1[163]*Gx2[135] + Gx1[164]*Gx2[149] + Gx1[165]*Gx2[163] + Gx1[166]*Gx2[177] + Gx1[167]*Gx2[191];
Gx3[164] = + Gx1[154]*Gx2[10] + Gx1[155]*Gx2[24] + Gx1[156]*Gx2[38] + Gx1[157]*Gx2[52] + Gx1[158]*Gx2[66] + Gx1[159]*Gx2[80] + Gx1[160]*Gx2[94] + Gx1[161]*Gx2[108] + Gx1[162]*Gx2[122] + Gx1[163]*Gx2[136] + Gx1[164]*Gx2[150] + Gx1[165]*Gx2[164] + Gx1[166]*Gx2[178] + Gx1[167]*Gx2[192];
Gx3[165] = + Gx1[154]*Gx2[11] + Gx1[155]*Gx2[25] + Gx1[156]*Gx2[39] + Gx1[157]*Gx2[53] + Gx1[158]*Gx2[67] + Gx1[159]*Gx2[81] + Gx1[160]*Gx2[95] + Gx1[161]*Gx2[109] + Gx1[162]*Gx2[123] + Gx1[163]*Gx2[137] + Gx1[164]*Gx2[151] + Gx1[165]*Gx2[165] + Gx1[166]*Gx2[179] + Gx1[167]*Gx2[193];
Gx3[166] = + Gx1[154]*Gx2[12] + Gx1[155]*Gx2[26] + Gx1[156]*Gx2[40] + Gx1[157]*Gx2[54] + Gx1[158]*Gx2[68] + Gx1[159]*Gx2[82] + Gx1[160]*Gx2[96] + Gx1[161]*Gx2[110] + Gx1[162]*Gx2[124] + Gx1[163]*Gx2[138] + Gx1[164]*Gx2[152] + Gx1[165]*Gx2[166] + Gx1[166]*Gx2[180] + Gx1[167]*Gx2[194];
Gx3[167] = + Gx1[154]*Gx2[13] + Gx1[155]*Gx2[27] + Gx1[156]*Gx2[41] + Gx1[157]*Gx2[55] + Gx1[158]*Gx2[69] + Gx1[159]*Gx2[83] + Gx1[160]*Gx2[97] + Gx1[161]*Gx2[111] + Gx1[162]*Gx2[125] + Gx1[163]*Gx2[139] + Gx1[164]*Gx2[153] + Gx1[165]*Gx2[167] + Gx1[166]*Gx2[181] + Gx1[167]*Gx2[195];
Gx3[168] = + Gx1[168]*Gx2[0] + Gx1[169]*Gx2[14] + Gx1[170]*Gx2[28] + Gx1[171]*Gx2[42] + Gx1[172]*Gx2[56] + Gx1[173]*Gx2[70] + Gx1[174]*Gx2[84] + Gx1[175]*Gx2[98] + Gx1[176]*Gx2[112] + Gx1[177]*Gx2[126] + Gx1[178]*Gx2[140] + Gx1[179]*Gx2[154] + Gx1[180]*Gx2[168] + Gx1[181]*Gx2[182];
Gx3[169] = + Gx1[168]*Gx2[1] + Gx1[169]*Gx2[15] + Gx1[170]*Gx2[29] + Gx1[171]*Gx2[43] + Gx1[172]*Gx2[57] + Gx1[173]*Gx2[71] + Gx1[174]*Gx2[85] + Gx1[175]*Gx2[99] + Gx1[176]*Gx2[113] + Gx1[177]*Gx2[127] + Gx1[178]*Gx2[141] + Gx1[179]*Gx2[155] + Gx1[180]*Gx2[169] + Gx1[181]*Gx2[183];
Gx3[170] = + Gx1[168]*Gx2[2] + Gx1[169]*Gx2[16] + Gx1[170]*Gx2[30] + Gx1[171]*Gx2[44] + Gx1[172]*Gx2[58] + Gx1[173]*Gx2[72] + Gx1[174]*Gx2[86] + Gx1[175]*Gx2[100] + Gx1[176]*Gx2[114] + Gx1[177]*Gx2[128] + Gx1[178]*Gx2[142] + Gx1[179]*Gx2[156] + Gx1[180]*Gx2[170] + Gx1[181]*Gx2[184];
Gx3[171] = + Gx1[168]*Gx2[3] + Gx1[169]*Gx2[17] + Gx1[170]*Gx2[31] + Gx1[171]*Gx2[45] + Gx1[172]*Gx2[59] + Gx1[173]*Gx2[73] + Gx1[174]*Gx2[87] + Gx1[175]*Gx2[101] + Gx1[176]*Gx2[115] + Gx1[177]*Gx2[129] + Gx1[178]*Gx2[143] + Gx1[179]*Gx2[157] + Gx1[180]*Gx2[171] + Gx1[181]*Gx2[185];
Gx3[172] = + Gx1[168]*Gx2[4] + Gx1[169]*Gx2[18] + Gx1[170]*Gx2[32] + Gx1[171]*Gx2[46] + Gx1[172]*Gx2[60] + Gx1[173]*Gx2[74] + Gx1[174]*Gx2[88] + Gx1[175]*Gx2[102] + Gx1[176]*Gx2[116] + Gx1[177]*Gx2[130] + Gx1[178]*Gx2[144] + Gx1[179]*Gx2[158] + Gx1[180]*Gx2[172] + Gx1[181]*Gx2[186];
Gx3[173] = + Gx1[168]*Gx2[5] + Gx1[169]*Gx2[19] + Gx1[170]*Gx2[33] + Gx1[171]*Gx2[47] + Gx1[172]*Gx2[61] + Gx1[173]*Gx2[75] + Gx1[174]*Gx2[89] + Gx1[175]*Gx2[103] + Gx1[176]*Gx2[117] + Gx1[177]*Gx2[131] + Gx1[178]*Gx2[145] + Gx1[179]*Gx2[159] + Gx1[180]*Gx2[173] + Gx1[181]*Gx2[187];
Gx3[174] = + Gx1[168]*Gx2[6] + Gx1[169]*Gx2[20] + Gx1[170]*Gx2[34] + Gx1[171]*Gx2[48] + Gx1[172]*Gx2[62] + Gx1[173]*Gx2[76] + Gx1[174]*Gx2[90] + Gx1[175]*Gx2[104] + Gx1[176]*Gx2[118] + Gx1[177]*Gx2[132] + Gx1[178]*Gx2[146] + Gx1[179]*Gx2[160] + Gx1[180]*Gx2[174] + Gx1[181]*Gx2[188];
Gx3[175] = + Gx1[168]*Gx2[7] + Gx1[169]*Gx2[21] + Gx1[170]*Gx2[35] + Gx1[171]*Gx2[49] + Gx1[172]*Gx2[63] + Gx1[173]*Gx2[77] + Gx1[174]*Gx2[91] + Gx1[175]*Gx2[105] + Gx1[176]*Gx2[119] + Gx1[177]*Gx2[133] + Gx1[178]*Gx2[147] + Gx1[179]*Gx2[161] + Gx1[180]*Gx2[175] + Gx1[181]*Gx2[189];
Gx3[176] = + Gx1[168]*Gx2[8] + Gx1[169]*Gx2[22] + Gx1[170]*Gx2[36] + Gx1[171]*Gx2[50] + Gx1[172]*Gx2[64] + Gx1[173]*Gx2[78] + Gx1[174]*Gx2[92] + Gx1[175]*Gx2[106] + Gx1[176]*Gx2[120] + Gx1[177]*Gx2[134] + Gx1[178]*Gx2[148] + Gx1[179]*Gx2[162] + Gx1[180]*Gx2[176] + Gx1[181]*Gx2[190];
Gx3[177] = + Gx1[168]*Gx2[9] + Gx1[169]*Gx2[23] + Gx1[170]*Gx2[37] + Gx1[171]*Gx2[51] + Gx1[172]*Gx2[65] + Gx1[173]*Gx2[79] + Gx1[174]*Gx2[93] + Gx1[175]*Gx2[107] + Gx1[176]*Gx2[121] + Gx1[177]*Gx2[135] + Gx1[178]*Gx2[149] + Gx1[179]*Gx2[163] + Gx1[180]*Gx2[177] + Gx1[181]*Gx2[191];
Gx3[178] = + Gx1[168]*Gx2[10] + Gx1[169]*Gx2[24] + Gx1[170]*Gx2[38] + Gx1[171]*Gx2[52] + Gx1[172]*Gx2[66] + Gx1[173]*Gx2[80] + Gx1[174]*Gx2[94] + Gx1[175]*Gx2[108] + Gx1[176]*Gx2[122] + Gx1[177]*Gx2[136] + Gx1[178]*Gx2[150] + Gx1[179]*Gx2[164] + Gx1[180]*Gx2[178] + Gx1[181]*Gx2[192];
Gx3[179] = + Gx1[168]*Gx2[11] + Gx1[169]*Gx2[25] + Gx1[170]*Gx2[39] + Gx1[171]*Gx2[53] + Gx1[172]*Gx2[67] + Gx1[173]*Gx2[81] + Gx1[174]*Gx2[95] + Gx1[175]*Gx2[109] + Gx1[176]*Gx2[123] + Gx1[177]*Gx2[137] + Gx1[178]*Gx2[151] + Gx1[179]*Gx2[165] + Gx1[180]*Gx2[179] + Gx1[181]*Gx2[193];
Gx3[180] = + Gx1[168]*Gx2[12] + Gx1[169]*Gx2[26] + Gx1[170]*Gx2[40] + Gx1[171]*Gx2[54] + Gx1[172]*Gx2[68] + Gx1[173]*Gx2[82] + Gx1[174]*Gx2[96] + Gx1[175]*Gx2[110] + Gx1[176]*Gx2[124] + Gx1[177]*Gx2[138] + Gx1[178]*Gx2[152] + Gx1[179]*Gx2[166] + Gx1[180]*Gx2[180] + Gx1[181]*Gx2[194];
Gx3[181] = + Gx1[168]*Gx2[13] + Gx1[169]*Gx2[27] + Gx1[170]*Gx2[41] + Gx1[171]*Gx2[55] + Gx1[172]*Gx2[69] + Gx1[173]*Gx2[83] + Gx1[174]*Gx2[97] + Gx1[175]*Gx2[111] + Gx1[176]*Gx2[125] + Gx1[177]*Gx2[139] + Gx1[178]*Gx2[153] + Gx1[179]*Gx2[167] + Gx1[180]*Gx2[181] + Gx1[181]*Gx2[195];
Gx3[182] = + Gx1[182]*Gx2[0] + Gx1[183]*Gx2[14] + Gx1[184]*Gx2[28] + Gx1[185]*Gx2[42] + Gx1[186]*Gx2[56] + Gx1[187]*Gx2[70] + Gx1[188]*Gx2[84] + Gx1[189]*Gx2[98] + Gx1[190]*Gx2[112] + Gx1[191]*Gx2[126] + Gx1[192]*Gx2[140] + Gx1[193]*Gx2[154] + Gx1[194]*Gx2[168] + Gx1[195]*Gx2[182];
Gx3[183] = + Gx1[182]*Gx2[1] + Gx1[183]*Gx2[15] + Gx1[184]*Gx2[29] + Gx1[185]*Gx2[43] + Gx1[186]*Gx2[57] + Gx1[187]*Gx2[71] + Gx1[188]*Gx2[85] + Gx1[189]*Gx2[99] + Gx1[190]*Gx2[113] + Gx1[191]*Gx2[127] + Gx1[192]*Gx2[141] + Gx1[193]*Gx2[155] + Gx1[194]*Gx2[169] + Gx1[195]*Gx2[183];
Gx3[184] = + Gx1[182]*Gx2[2] + Gx1[183]*Gx2[16] + Gx1[184]*Gx2[30] + Gx1[185]*Gx2[44] + Gx1[186]*Gx2[58] + Gx1[187]*Gx2[72] + Gx1[188]*Gx2[86] + Gx1[189]*Gx2[100] + Gx1[190]*Gx2[114] + Gx1[191]*Gx2[128] + Gx1[192]*Gx2[142] + Gx1[193]*Gx2[156] + Gx1[194]*Gx2[170] + Gx1[195]*Gx2[184];
Gx3[185] = + Gx1[182]*Gx2[3] + Gx1[183]*Gx2[17] + Gx1[184]*Gx2[31] + Gx1[185]*Gx2[45] + Gx1[186]*Gx2[59] + Gx1[187]*Gx2[73] + Gx1[188]*Gx2[87] + Gx1[189]*Gx2[101] + Gx1[190]*Gx2[115] + Gx1[191]*Gx2[129] + Gx1[192]*Gx2[143] + Gx1[193]*Gx2[157] + Gx1[194]*Gx2[171] + Gx1[195]*Gx2[185];
Gx3[186] = + Gx1[182]*Gx2[4] + Gx1[183]*Gx2[18] + Gx1[184]*Gx2[32] + Gx1[185]*Gx2[46] + Gx1[186]*Gx2[60] + Gx1[187]*Gx2[74] + Gx1[188]*Gx2[88] + Gx1[189]*Gx2[102] + Gx1[190]*Gx2[116] + Gx1[191]*Gx2[130] + Gx1[192]*Gx2[144] + Gx1[193]*Gx2[158] + Gx1[194]*Gx2[172] + Gx1[195]*Gx2[186];
Gx3[187] = + Gx1[182]*Gx2[5] + Gx1[183]*Gx2[19] + Gx1[184]*Gx2[33] + Gx1[185]*Gx2[47] + Gx1[186]*Gx2[61] + Gx1[187]*Gx2[75] + Gx1[188]*Gx2[89] + Gx1[189]*Gx2[103] + Gx1[190]*Gx2[117] + Gx1[191]*Gx2[131] + Gx1[192]*Gx2[145] + Gx1[193]*Gx2[159] + Gx1[194]*Gx2[173] + Gx1[195]*Gx2[187];
Gx3[188] = + Gx1[182]*Gx2[6] + Gx1[183]*Gx2[20] + Gx1[184]*Gx2[34] + Gx1[185]*Gx2[48] + Gx1[186]*Gx2[62] + Gx1[187]*Gx2[76] + Gx1[188]*Gx2[90] + Gx1[189]*Gx2[104] + Gx1[190]*Gx2[118] + Gx1[191]*Gx2[132] + Gx1[192]*Gx2[146] + Gx1[193]*Gx2[160] + Gx1[194]*Gx2[174] + Gx1[195]*Gx2[188];
Gx3[189] = + Gx1[182]*Gx2[7] + Gx1[183]*Gx2[21] + Gx1[184]*Gx2[35] + Gx1[185]*Gx2[49] + Gx1[186]*Gx2[63] + Gx1[187]*Gx2[77] + Gx1[188]*Gx2[91] + Gx1[189]*Gx2[105] + Gx1[190]*Gx2[119] + Gx1[191]*Gx2[133] + Gx1[192]*Gx2[147] + Gx1[193]*Gx2[161] + Gx1[194]*Gx2[175] + Gx1[195]*Gx2[189];
Gx3[190] = + Gx1[182]*Gx2[8] + Gx1[183]*Gx2[22] + Gx1[184]*Gx2[36] + Gx1[185]*Gx2[50] + Gx1[186]*Gx2[64] + Gx1[187]*Gx2[78] + Gx1[188]*Gx2[92] + Gx1[189]*Gx2[106] + Gx1[190]*Gx2[120] + Gx1[191]*Gx2[134] + Gx1[192]*Gx2[148] + Gx1[193]*Gx2[162] + Gx1[194]*Gx2[176] + Gx1[195]*Gx2[190];
Gx3[191] = + Gx1[182]*Gx2[9] + Gx1[183]*Gx2[23] + Gx1[184]*Gx2[37] + Gx1[185]*Gx2[51] + Gx1[186]*Gx2[65] + Gx1[187]*Gx2[79] + Gx1[188]*Gx2[93] + Gx1[189]*Gx2[107] + Gx1[190]*Gx2[121] + Gx1[191]*Gx2[135] + Gx1[192]*Gx2[149] + Gx1[193]*Gx2[163] + Gx1[194]*Gx2[177] + Gx1[195]*Gx2[191];
Gx3[192] = + Gx1[182]*Gx2[10] + Gx1[183]*Gx2[24] + Gx1[184]*Gx2[38] + Gx1[185]*Gx2[52] + Gx1[186]*Gx2[66] + Gx1[187]*Gx2[80] + Gx1[188]*Gx2[94] + Gx1[189]*Gx2[108] + Gx1[190]*Gx2[122] + Gx1[191]*Gx2[136] + Gx1[192]*Gx2[150] + Gx1[193]*Gx2[164] + Gx1[194]*Gx2[178] + Gx1[195]*Gx2[192];
Gx3[193] = + Gx1[182]*Gx2[11] + Gx1[183]*Gx2[25] + Gx1[184]*Gx2[39] + Gx1[185]*Gx2[53] + Gx1[186]*Gx2[67] + Gx1[187]*Gx2[81] + Gx1[188]*Gx2[95] + Gx1[189]*Gx2[109] + Gx1[190]*Gx2[123] + Gx1[191]*Gx2[137] + Gx1[192]*Gx2[151] + Gx1[193]*Gx2[165] + Gx1[194]*Gx2[179] + Gx1[195]*Gx2[193];
Gx3[194] = + Gx1[182]*Gx2[12] + Gx1[183]*Gx2[26] + Gx1[184]*Gx2[40] + Gx1[185]*Gx2[54] + Gx1[186]*Gx2[68] + Gx1[187]*Gx2[82] + Gx1[188]*Gx2[96] + Gx1[189]*Gx2[110] + Gx1[190]*Gx2[124] + Gx1[191]*Gx2[138] + Gx1[192]*Gx2[152] + Gx1[193]*Gx2[166] + Gx1[194]*Gx2[180] + Gx1[195]*Gx2[194];
Gx3[195] = + Gx1[182]*Gx2[13] + Gx1[183]*Gx2[27] + Gx1[184]*Gx2[41] + Gx1[185]*Gx2[55] + Gx1[186]*Gx2[69] + Gx1[187]*Gx2[83] + Gx1[188]*Gx2[97] + Gx1[189]*Gx2[111] + Gx1[190]*Gx2[125] + Gx1[191]*Gx2[139] + Gx1[192]*Gx2[153] + Gx1[193]*Gx2[167] + Gx1[194]*Gx2[181] + Gx1[195]*Gx2[195];
}

void nmpc_multGxGu( real_t* const Gx1, real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = + Gx1[0]*Gu1[0] + Gx1[1]*Gu1[4] + Gx1[2]*Gu1[8] + Gx1[3]*Gu1[12] + Gx1[4]*Gu1[16] + Gx1[5]*Gu1[20] + Gx1[6]*Gu1[24] + Gx1[7]*Gu1[28] + Gx1[8]*Gu1[32] + Gx1[9]*Gu1[36] + Gx1[10]*Gu1[40] + Gx1[11]*Gu1[44] + Gx1[12]*Gu1[48] + Gx1[13]*Gu1[52];
Gu2[1] = + Gx1[0]*Gu1[1] + Gx1[1]*Gu1[5] + Gx1[2]*Gu1[9] + Gx1[3]*Gu1[13] + Gx1[4]*Gu1[17] + Gx1[5]*Gu1[21] + Gx1[6]*Gu1[25] + Gx1[7]*Gu1[29] + Gx1[8]*Gu1[33] + Gx1[9]*Gu1[37] + Gx1[10]*Gu1[41] + Gx1[11]*Gu1[45] + Gx1[12]*Gu1[49] + Gx1[13]*Gu1[53];
Gu2[2] = + Gx1[0]*Gu1[2] + Gx1[1]*Gu1[6] + Gx1[2]*Gu1[10] + Gx1[3]*Gu1[14] + Gx1[4]*Gu1[18] + Gx1[5]*Gu1[22] + Gx1[6]*Gu1[26] + Gx1[7]*Gu1[30] + Gx1[8]*Gu1[34] + Gx1[9]*Gu1[38] + Gx1[10]*Gu1[42] + Gx1[11]*Gu1[46] + Gx1[12]*Gu1[50] + Gx1[13]*Gu1[54];
Gu2[3] = + Gx1[0]*Gu1[3] + Gx1[1]*Gu1[7] + Gx1[2]*Gu1[11] + Gx1[3]*Gu1[15] + Gx1[4]*Gu1[19] + Gx1[5]*Gu1[23] + Gx1[6]*Gu1[27] + Gx1[7]*Gu1[31] + Gx1[8]*Gu1[35] + Gx1[9]*Gu1[39] + Gx1[10]*Gu1[43] + Gx1[11]*Gu1[47] + Gx1[12]*Gu1[51] + Gx1[13]*Gu1[55];
Gu2[4] = + Gx1[14]*Gu1[0] + Gx1[15]*Gu1[4] + Gx1[16]*Gu1[8] + Gx1[17]*Gu1[12] + Gx1[18]*Gu1[16] + Gx1[19]*Gu1[20] + Gx1[20]*Gu1[24] + Gx1[21]*Gu1[28] + Gx1[22]*Gu1[32] + Gx1[23]*Gu1[36] + Gx1[24]*Gu1[40] + Gx1[25]*Gu1[44] + Gx1[26]*Gu1[48] + Gx1[27]*Gu1[52];
Gu2[5] = + Gx1[14]*Gu1[1] + Gx1[15]*Gu1[5] + Gx1[16]*Gu1[9] + Gx1[17]*Gu1[13] + Gx1[18]*Gu1[17] + Gx1[19]*Gu1[21] + Gx1[20]*Gu1[25] + Gx1[21]*Gu1[29] + Gx1[22]*Gu1[33] + Gx1[23]*Gu1[37] + Gx1[24]*Gu1[41] + Gx1[25]*Gu1[45] + Gx1[26]*Gu1[49] + Gx1[27]*Gu1[53];
Gu2[6] = + Gx1[14]*Gu1[2] + Gx1[15]*Gu1[6] + Gx1[16]*Gu1[10] + Gx1[17]*Gu1[14] + Gx1[18]*Gu1[18] + Gx1[19]*Gu1[22] + Gx1[20]*Gu1[26] + Gx1[21]*Gu1[30] + Gx1[22]*Gu1[34] + Gx1[23]*Gu1[38] + Gx1[24]*Gu1[42] + Gx1[25]*Gu1[46] + Gx1[26]*Gu1[50] + Gx1[27]*Gu1[54];
Gu2[7] = + Gx1[14]*Gu1[3] + Gx1[15]*Gu1[7] + Gx1[16]*Gu1[11] + Gx1[17]*Gu1[15] + Gx1[18]*Gu1[19] + Gx1[19]*Gu1[23] + Gx1[20]*Gu1[27] + Gx1[21]*Gu1[31] + Gx1[22]*Gu1[35] + Gx1[23]*Gu1[39] + Gx1[24]*Gu1[43] + Gx1[25]*Gu1[47] + Gx1[26]*Gu1[51] + Gx1[27]*Gu1[55];
Gu2[8] = + Gx1[28]*Gu1[0] + Gx1[29]*Gu1[4] + Gx1[30]*Gu1[8] + Gx1[31]*Gu1[12] + Gx1[32]*Gu1[16] + Gx1[33]*Gu1[20] + Gx1[34]*Gu1[24] + Gx1[35]*Gu1[28] + Gx1[36]*Gu1[32] + Gx1[37]*Gu1[36] + Gx1[38]*Gu1[40] + Gx1[39]*Gu1[44] + Gx1[40]*Gu1[48] + Gx1[41]*Gu1[52];
Gu2[9] = + Gx1[28]*Gu1[1] + Gx1[29]*Gu1[5] + Gx1[30]*Gu1[9] + Gx1[31]*Gu1[13] + Gx1[32]*Gu1[17] + Gx1[33]*Gu1[21] + Gx1[34]*Gu1[25] + Gx1[35]*Gu1[29] + Gx1[36]*Gu1[33] + Gx1[37]*Gu1[37] + Gx1[38]*Gu1[41] + Gx1[39]*Gu1[45] + Gx1[40]*Gu1[49] + Gx1[41]*Gu1[53];
Gu2[10] = + Gx1[28]*Gu1[2] + Gx1[29]*Gu1[6] + Gx1[30]*Gu1[10] + Gx1[31]*Gu1[14] + Gx1[32]*Gu1[18] + Gx1[33]*Gu1[22] + Gx1[34]*Gu1[26] + Gx1[35]*Gu1[30] + Gx1[36]*Gu1[34] + Gx1[37]*Gu1[38] + Gx1[38]*Gu1[42] + Gx1[39]*Gu1[46] + Gx1[40]*Gu1[50] + Gx1[41]*Gu1[54];
Gu2[11] = + Gx1[28]*Gu1[3] + Gx1[29]*Gu1[7] + Gx1[30]*Gu1[11] + Gx1[31]*Gu1[15] + Gx1[32]*Gu1[19] + Gx1[33]*Gu1[23] + Gx1[34]*Gu1[27] + Gx1[35]*Gu1[31] + Gx1[36]*Gu1[35] + Gx1[37]*Gu1[39] + Gx1[38]*Gu1[43] + Gx1[39]*Gu1[47] + Gx1[40]*Gu1[51] + Gx1[41]*Gu1[55];
Gu2[12] = + Gx1[42]*Gu1[0] + Gx1[43]*Gu1[4] + Gx1[44]*Gu1[8] + Gx1[45]*Gu1[12] + Gx1[46]*Gu1[16] + Gx1[47]*Gu1[20] + Gx1[48]*Gu1[24] + Gx1[49]*Gu1[28] + Gx1[50]*Gu1[32] + Gx1[51]*Gu1[36] + Gx1[52]*Gu1[40] + Gx1[53]*Gu1[44] + Gx1[54]*Gu1[48] + Gx1[55]*Gu1[52];
Gu2[13] = + Gx1[42]*Gu1[1] + Gx1[43]*Gu1[5] + Gx1[44]*Gu1[9] + Gx1[45]*Gu1[13] + Gx1[46]*Gu1[17] + Gx1[47]*Gu1[21] + Gx1[48]*Gu1[25] + Gx1[49]*Gu1[29] + Gx1[50]*Gu1[33] + Gx1[51]*Gu1[37] + Gx1[52]*Gu1[41] + Gx1[53]*Gu1[45] + Gx1[54]*Gu1[49] + Gx1[55]*Gu1[53];
Gu2[14] = + Gx1[42]*Gu1[2] + Gx1[43]*Gu1[6] + Gx1[44]*Gu1[10] + Gx1[45]*Gu1[14] + Gx1[46]*Gu1[18] + Gx1[47]*Gu1[22] + Gx1[48]*Gu1[26] + Gx1[49]*Gu1[30] + Gx1[50]*Gu1[34] + Gx1[51]*Gu1[38] + Gx1[52]*Gu1[42] + Gx1[53]*Gu1[46] + Gx1[54]*Gu1[50] + Gx1[55]*Gu1[54];
Gu2[15] = + Gx1[42]*Gu1[3] + Gx1[43]*Gu1[7] + Gx1[44]*Gu1[11] + Gx1[45]*Gu1[15] + Gx1[46]*Gu1[19] + Gx1[47]*Gu1[23] + Gx1[48]*Gu1[27] + Gx1[49]*Gu1[31] + Gx1[50]*Gu1[35] + Gx1[51]*Gu1[39] + Gx1[52]*Gu1[43] + Gx1[53]*Gu1[47] + Gx1[54]*Gu1[51] + Gx1[55]*Gu1[55];
Gu2[16] = + Gx1[56]*Gu1[0] + Gx1[57]*Gu1[4] + Gx1[58]*Gu1[8] + Gx1[59]*Gu1[12] + Gx1[60]*Gu1[16] + Gx1[61]*Gu1[20] + Gx1[62]*Gu1[24] + Gx1[63]*Gu1[28] + Gx1[64]*Gu1[32] + Gx1[65]*Gu1[36] + Gx1[66]*Gu1[40] + Gx1[67]*Gu1[44] + Gx1[68]*Gu1[48] + Gx1[69]*Gu1[52];
Gu2[17] = + Gx1[56]*Gu1[1] + Gx1[57]*Gu1[5] + Gx1[58]*Gu1[9] + Gx1[59]*Gu1[13] + Gx1[60]*Gu1[17] + Gx1[61]*Gu1[21] + Gx1[62]*Gu1[25] + Gx1[63]*Gu1[29] + Gx1[64]*Gu1[33] + Gx1[65]*Gu1[37] + Gx1[66]*Gu1[41] + Gx1[67]*Gu1[45] + Gx1[68]*Gu1[49] + Gx1[69]*Gu1[53];
Gu2[18] = + Gx1[56]*Gu1[2] + Gx1[57]*Gu1[6] + Gx1[58]*Gu1[10] + Gx1[59]*Gu1[14] + Gx1[60]*Gu1[18] + Gx1[61]*Gu1[22] + Gx1[62]*Gu1[26] + Gx1[63]*Gu1[30] + Gx1[64]*Gu1[34] + Gx1[65]*Gu1[38] + Gx1[66]*Gu1[42] + Gx1[67]*Gu1[46] + Gx1[68]*Gu1[50] + Gx1[69]*Gu1[54];
Gu2[19] = + Gx1[56]*Gu1[3] + Gx1[57]*Gu1[7] + Gx1[58]*Gu1[11] + Gx1[59]*Gu1[15] + Gx1[60]*Gu1[19] + Gx1[61]*Gu1[23] + Gx1[62]*Gu1[27] + Gx1[63]*Gu1[31] + Gx1[64]*Gu1[35] + Gx1[65]*Gu1[39] + Gx1[66]*Gu1[43] + Gx1[67]*Gu1[47] + Gx1[68]*Gu1[51] + Gx1[69]*Gu1[55];
Gu2[20] = + Gx1[70]*Gu1[0] + Gx1[71]*Gu1[4] + Gx1[72]*Gu1[8] + Gx1[73]*Gu1[12] + Gx1[74]*Gu1[16] + Gx1[75]*Gu1[20] + Gx1[76]*Gu1[24] + Gx1[77]*Gu1[28] + Gx1[78]*Gu1[32] + Gx1[79]*Gu1[36] + Gx1[80]*Gu1[40] + Gx1[81]*Gu1[44] + Gx1[82]*Gu1[48] + Gx1[83]*Gu1[52];
Gu2[21] = + Gx1[70]*Gu1[1] + Gx1[71]*Gu1[5] + Gx1[72]*Gu1[9] + Gx1[73]*Gu1[13] + Gx1[74]*Gu1[17] + Gx1[75]*Gu1[21] + Gx1[76]*Gu1[25] + Gx1[77]*Gu1[29] + Gx1[78]*Gu1[33] + Gx1[79]*Gu1[37] + Gx1[80]*Gu1[41] + Gx1[81]*Gu1[45] + Gx1[82]*Gu1[49] + Gx1[83]*Gu1[53];
Gu2[22] = + Gx1[70]*Gu1[2] + Gx1[71]*Gu1[6] + Gx1[72]*Gu1[10] + Gx1[73]*Gu1[14] + Gx1[74]*Gu1[18] + Gx1[75]*Gu1[22] + Gx1[76]*Gu1[26] + Gx1[77]*Gu1[30] + Gx1[78]*Gu1[34] + Gx1[79]*Gu1[38] + Gx1[80]*Gu1[42] + Gx1[81]*Gu1[46] + Gx1[82]*Gu1[50] + Gx1[83]*Gu1[54];
Gu2[23] = + Gx1[70]*Gu1[3] + Gx1[71]*Gu1[7] + Gx1[72]*Gu1[11] + Gx1[73]*Gu1[15] + Gx1[74]*Gu1[19] + Gx1[75]*Gu1[23] + Gx1[76]*Gu1[27] + Gx1[77]*Gu1[31] + Gx1[78]*Gu1[35] + Gx1[79]*Gu1[39] + Gx1[80]*Gu1[43] + Gx1[81]*Gu1[47] + Gx1[82]*Gu1[51] + Gx1[83]*Gu1[55];
Gu2[24] = + Gx1[84]*Gu1[0] + Gx1[85]*Gu1[4] + Gx1[86]*Gu1[8] + Gx1[87]*Gu1[12] + Gx1[88]*Gu1[16] + Gx1[89]*Gu1[20] + Gx1[90]*Gu1[24] + Gx1[91]*Gu1[28] + Gx1[92]*Gu1[32] + Gx1[93]*Gu1[36] + Gx1[94]*Gu1[40] + Gx1[95]*Gu1[44] + Gx1[96]*Gu1[48] + Gx1[97]*Gu1[52];
Gu2[25] = + Gx1[84]*Gu1[1] + Gx1[85]*Gu1[5] + Gx1[86]*Gu1[9] + Gx1[87]*Gu1[13] + Gx1[88]*Gu1[17] + Gx1[89]*Gu1[21] + Gx1[90]*Gu1[25] + Gx1[91]*Gu1[29] + Gx1[92]*Gu1[33] + Gx1[93]*Gu1[37] + Gx1[94]*Gu1[41] + Gx1[95]*Gu1[45] + Gx1[96]*Gu1[49] + Gx1[97]*Gu1[53];
Gu2[26] = + Gx1[84]*Gu1[2] + Gx1[85]*Gu1[6] + Gx1[86]*Gu1[10] + Gx1[87]*Gu1[14] + Gx1[88]*Gu1[18] + Gx1[89]*Gu1[22] + Gx1[90]*Gu1[26] + Gx1[91]*Gu1[30] + Gx1[92]*Gu1[34] + Gx1[93]*Gu1[38] + Gx1[94]*Gu1[42] + Gx1[95]*Gu1[46] + Gx1[96]*Gu1[50] + Gx1[97]*Gu1[54];
Gu2[27] = + Gx1[84]*Gu1[3] + Gx1[85]*Gu1[7] + Gx1[86]*Gu1[11] + Gx1[87]*Gu1[15] + Gx1[88]*Gu1[19] + Gx1[89]*Gu1[23] + Gx1[90]*Gu1[27] + Gx1[91]*Gu1[31] + Gx1[92]*Gu1[35] + Gx1[93]*Gu1[39] + Gx1[94]*Gu1[43] + Gx1[95]*Gu1[47] + Gx1[96]*Gu1[51] + Gx1[97]*Gu1[55];
Gu2[28] = + Gx1[98]*Gu1[0] + Gx1[99]*Gu1[4] + Gx1[100]*Gu1[8] + Gx1[101]*Gu1[12] + Gx1[102]*Gu1[16] + Gx1[103]*Gu1[20] + Gx1[104]*Gu1[24] + Gx1[105]*Gu1[28] + Gx1[106]*Gu1[32] + Gx1[107]*Gu1[36] + Gx1[108]*Gu1[40] + Gx1[109]*Gu1[44] + Gx1[110]*Gu1[48] + Gx1[111]*Gu1[52];
Gu2[29] = + Gx1[98]*Gu1[1] + Gx1[99]*Gu1[5] + Gx1[100]*Gu1[9] + Gx1[101]*Gu1[13] + Gx1[102]*Gu1[17] + Gx1[103]*Gu1[21] + Gx1[104]*Gu1[25] + Gx1[105]*Gu1[29] + Gx1[106]*Gu1[33] + Gx1[107]*Gu1[37] + Gx1[108]*Gu1[41] + Gx1[109]*Gu1[45] + Gx1[110]*Gu1[49] + Gx1[111]*Gu1[53];
Gu2[30] = + Gx1[98]*Gu1[2] + Gx1[99]*Gu1[6] + Gx1[100]*Gu1[10] + Gx1[101]*Gu1[14] + Gx1[102]*Gu1[18] + Gx1[103]*Gu1[22] + Gx1[104]*Gu1[26] + Gx1[105]*Gu1[30] + Gx1[106]*Gu1[34] + Gx1[107]*Gu1[38] + Gx1[108]*Gu1[42] + Gx1[109]*Gu1[46] + Gx1[110]*Gu1[50] + Gx1[111]*Gu1[54];
Gu2[31] = + Gx1[98]*Gu1[3] + Gx1[99]*Gu1[7] + Gx1[100]*Gu1[11] + Gx1[101]*Gu1[15] + Gx1[102]*Gu1[19] + Gx1[103]*Gu1[23] + Gx1[104]*Gu1[27] + Gx1[105]*Gu1[31] + Gx1[106]*Gu1[35] + Gx1[107]*Gu1[39] + Gx1[108]*Gu1[43] + Gx1[109]*Gu1[47] + Gx1[110]*Gu1[51] + Gx1[111]*Gu1[55];
Gu2[32] = + Gx1[112]*Gu1[0] + Gx1[113]*Gu1[4] + Gx1[114]*Gu1[8] + Gx1[115]*Gu1[12] + Gx1[116]*Gu1[16] + Gx1[117]*Gu1[20] + Gx1[118]*Gu1[24] + Gx1[119]*Gu1[28] + Gx1[120]*Gu1[32] + Gx1[121]*Gu1[36] + Gx1[122]*Gu1[40] + Gx1[123]*Gu1[44] + Gx1[124]*Gu1[48] + Gx1[125]*Gu1[52];
Gu2[33] = + Gx1[112]*Gu1[1] + Gx1[113]*Gu1[5] + Gx1[114]*Gu1[9] + Gx1[115]*Gu1[13] + Gx1[116]*Gu1[17] + Gx1[117]*Gu1[21] + Gx1[118]*Gu1[25] + Gx1[119]*Gu1[29] + Gx1[120]*Gu1[33] + Gx1[121]*Gu1[37] + Gx1[122]*Gu1[41] + Gx1[123]*Gu1[45] + Gx1[124]*Gu1[49] + Gx1[125]*Gu1[53];
Gu2[34] = + Gx1[112]*Gu1[2] + Gx1[113]*Gu1[6] + Gx1[114]*Gu1[10] + Gx1[115]*Gu1[14] + Gx1[116]*Gu1[18] + Gx1[117]*Gu1[22] + Gx1[118]*Gu1[26] + Gx1[119]*Gu1[30] + Gx1[120]*Gu1[34] + Gx1[121]*Gu1[38] + Gx1[122]*Gu1[42] + Gx1[123]*Gu1[46] + Gx1[124]*Gu1[50] + Gx1[125]*Gu1[54];
Gu2[35] = + Gx1[112]*Gu1[3] + Gx1[113]*Gu1[7] + Gx1[114]*Gu1[11] + Gx1[115]*Gu1[15] + Gx1[116]*Gu1[19] + Gx1[117]*Gu1[23] + Gx1[118]*Gu1[27] + Gx1[119]*Gu1[31] + Gx1[120]*Gu1[35] + Gx1[121]*Gu1[39] + Gx1[122]*Gu1[43] + Gx1[123]*Gu1[47] + Gx1[124]*Gu1[51] + Gx1[125]*Gu1[55];
Gu2[36] = + Gx1[126]*Gu1[0] + Gx1[127]*Gu1[4] + Gx1[128]*Gu1[8] + Gx1[129]*Gu1[12] + Gx1[130]*Gu1[16] + Gx1[131]*Gu1[20] + Gx1[132]*Gu1[24] + Gx1[133]*Gu1[28] + Gx1[134]*Gu1[32] + Gx1[135]*Gu1[36] + Gx1[136]*Gu1[40] + Gx1[137]*Gu1[44] + Gx1[138]*Gu1[48] + Gx1[139]*Gu1[52];
Gu2[37] = + Gx1[126]*Gu1[1] + Gx1[127]*Gu1[5] + Gx1[128]*Gu1[9] + Gx1[129]*Gu1[13] + Gx1[130]*Gu1[17] + Gx1[131]*Gu1[21] + Gx1[132]*Gu1[25] + Gx1[133]*Gu1[29] + Gx1[134]*Gu1[33] + Gx1[135]*Gu1[37] + Gx1[136]*Gu1[41] + Gx1[137]*Gu1[45] + Gx1[138]*Gu1[49] + Gx1[139]*Gu1[53];
Gu2[38] = + Gx1[126]*Gu1[2] + Gx1[127]*Gu1[6] + Gx1[128]*Gu1[10] + Gx1[129]*Gu1[14] + Gx1[130]*Gu1[18] + Gx1[131]*Gu1[22] + Gx1[132]*Gu1[26] + Gx1[133]*Gu1[30] + Gx1[134]*Gu1[34] + Gx1[135]*Gu1[38] + Gx1[136]*Gu1[42] + Gx1[137]*Gu1[46] + Gx1[138]*Gu1[50] + Gx1[139]*Gu1[54];
Gu2[39] = + Gx1[126]*Gu1[3] + Gx1[127]*Gu1[7] + Gx1[128]*Gu1[11] + Gx1[129]*Gu1[15] + Gx1[130]*Gu1[19] + Gx1[131]*Gu1[23] + Gx1[132]*Gu1[27] + Gx1[133]*Gu1[31] + Gx1[134]*Gu1[35] + Gx1[135]*Gu1[39] + Gx1[136]*Gu1[43] + Gx1[137]*Gu1[47] + Gx1[138]*Gu1[51] + Gx1[139]*Gu1[55];
Gu2[40] = + Gx1[140]*Gu1[0] + Gx1[141]*Gu1[4] + Gx1[142]*Gu1[8] + Gx1[143]*Gu1[12] + Gx1[144]*Gu1[16] + Gx1[145]*Gu1[20] + Gx1[146]*Gu1[24] + Gx1[147]*Gu1[28] + Gx1[148]*Gu1[32] + Gx1[149]*Gu1[36] + Gx1[150]*Gu1[40] + Gx1[151]*Gu1[44] + Gx1[152]*Gu1[48] + Gx1[153]*Gu1[52];
Gu2[41] = + Gx1[140]*Gu1[1] + Gx1[141]*Gu1[5] + Gx1[142]*Gu1[9] + Gx1[143]*Gu1[13] + Gx1[144]*Gu1[17] + Gx1[145]*Gu1[21] + Gx1[146]*Gu1[25] + Gx1[147]*Gu1[29] + Gx1[148]*Gu1[33] + Gx1[149]*Gu1[37] + Gx1[150]*Gu1[41] + Gx1[151]*Gu1[45] + Gx1[152]*Gu1[49] + Gx1[153]*Gu1[53];
Gu2[42] = + Gx1[140]*Gu1[2] + Gx1[141]*Gu1[6] + Gx1[142]*Gu1[10] + Gx1[143]*Gu1[14] + Gx1[144]*Gu1[18] + Gx1[145]*Gu1[22] + Gx1[146]*Gu1[26] + Gx1[147]*Gu1[30] + Gx1[148]*Gu1[34] + Gx1[149]*Gu1[38] + Gx1[150]*Gu1[42] + Gx1[151]*Gu1[46] + Gx1[152]*Gu1[50] + Gx1[153]*Gu1[54];
Gu2[43] = + Gx1[140]*Gu1[3] + Gx1[141]*Gu1[7] + Gx1[142]*Gu1[11] + Gx1[143]*Gu1[15] + Gx1[144]*Gu1[19] + Gx1[145]*Gu1[23] + Gx1[146]*Gu1[27] + Gx1[147]*Gu1[31] + Gx1[148]*Gu1[35] + Gx1[149]*Gu1[39] + Gx1[150]*Gu1[43] + Gx1[151]*Gu1[47] + Gx1[152]*Gu1[51] + Gx1[153]*Gu1[55];
Gu2[44] = + Gx1[154]*Gu1[0] + Gx1[155]*Gu1[4] + Gx1[156]*Gu1[8] + Gx1[157]*Gu1[12] + Gx1[158]*Gu1[16] + Gx1[159]*Gu1[20] + Gx1[160]*Gu1[24] + Gx1[161]*Gu1[28] + Gx1[162]*Gu1[32] + Gx1[163]*Gu1[36] + Gx1[164]*Gu1[40] + Gx1[165]*Gu1[44] + Gx1[166]*Gu1[48] + Gx1[167]*Gu1[52];
Gu2[45] = + Gx1[154]*Gu1[1] + Gx1[155]*Gu1[5] + Gx1[156]*Gu1[9] + Gx1[157]*Gu1[13] + Gx1[158]*Gu1[17] + Gx1[159]*Gu1[21] + Gx1[160]*Gu1[25] + Gx1[161]*Gu1[29] + Gx1[162]*Gu1[33] + Gx1[163]*Gu1[37] + Gx1[164]*Gu1[41] + Gx1[165]*Gu1[45] + Gx1[166]*Gu1[49] + Gx1[167]*Gu1[53];
Gu2[46] = + Gx1[154]*Gu1[2] + Gx1[155]*Gu1[6] + Gx1[156]*Gu1[10] + Gx1[157]*Gu1[14] + Gx1[158]*Gu1[18] + Gx1[159]*Gu1[22] + Gx1[160]*Gu1[26] + Gx1[161]*Gu1[30] + Gx1[162]*Gu1[34] + Gx1[163]*Gu1[38] + Gx1[164]*Gu1[42] + Gx1[165]*Gu1[46] + Gx1[166]*Gu1[50] + Gx1[167]*Gu1[54];
Gu2[47] = + Gx1[154]*Gu1[3] + Gx1[155]*Gu1[7] + Gx1[156]*Gu1[11] + Gx1[157]*Gu1[15] + Gx1[158]*Gu1[19] + Gx1[159]*Gu1[23] + Gx1[160]*Gu1[27] + Gx1[161]*Gu1[31] + Gx1[162]*Gu1[35] + Gx1[163]*Gu1[39] + Gx1[164]*Gu1[43] + Gx1[165]*Gu1[47] + Gx1[166]*Gu1[51] + Gx1[167]*Gu1[55];
Gu2[48] = + Gx1[168]*Gu1[0] + Gx1[169]*Gu1[4] + Gx1[170]*Gu1[8] + Gx1[171]*Gu1[12] + Gx1[172]*Gu1[16] + Gx1[173]*Gu1[20] + Gx1[174]*Gu1[24] + Gx1[175]*Gu1[28] + Gx1[176]*Gu1[32] + Gx1[177]*Gu1[36] + Gx1[178]*Gu1[40] + Gx1[179]*Gu1[44] + Gx1[180]*Gu1[48] + Gx1[181]*Gu1[52];
Gu2[49] = + Gx1[168]*Gu1[1] + Gx1[169]*Gu1[5] + Gx1[170]*Gu1[9] + Gx1[171]*Gu1[13] + Gx1[172]*Gu1[17] + Gx1[173]*Gu1[21] + Gx1[174]*Gu1[25] + Gx1[175]*Gu1[29] + Gx1[176]*Gu1[33] + Gx1[177]*Gu1[37] + Gx1[178]*Gu1[41] + Gx1[179]*Gu1[45] + Gx1[180]*Gu1[49] + Gx1[181]*Gu1[53];
Gu2[50] = + Gx1[168]*Gu1[2] + Gx1[169]*Gu1[6] + Gx1[170]*Gu1[10] + Gx1[171]*Gu1[14] + Gx1[172]*Gu1[18] + Gx1[173]*Gu1[22] + Gx1[174]*Gu1[26] + Gx1[175]*Gu1[30] + Gx1[176]*Gu1[34] + Gx1[177]*Gu1[38] + Gx1[178]*Gu1[42] + Gx1[179]*Gu1[46] + Gx1[180]*Gu1[50] + Gx1[181]*Gu1[54];
Gu2[51] = + Gx1[168]*Gu1[3] + Gx1[169]*Gu1[7] + Gx1[170]*Gu1[11] + Gx1[171]*Gu1[15] + Gx1[172]*Gu1[19] + Gx1[173]*Gu1[23] + Gx1[174]*Gu1[27] + Gx1[175]*Gu1[31] + Gx1[176]*Gu1[35] + Gx1[177]*Gu1[39] + Gx1[178]*Gu1[43] + Gx1[179]*Gu1[47] + Gx1[180]*Gu1[51] + Gx1[181]*Gu1[55];
Gu2[52] = + Gx1[182]*Gu1[0] + Gx1[183]*Gu1[4] + Gx1[184]*Gu1[8] + Gx1[185]*Gu1[12] + Gx1[186]*Gu1[16] + Gx1[187]*Gu1[20] + Gx1[188]*Gu1[24] + Gx1[189]*Gu1[28] + Gx1[190]*Gu1[32] + Gx1[191]*Gu1[36] + Gx1[192]*Gu1[40] + Gx1[193]*Gu1[44] + Gx1[194]*Gu1[48] + Gx1[195]*Gu1[52];
Gu2[53] = + Gx1[182]*Gu1[1] + Gx1[183]*Gu1[5] + Gx1[184]*Gu1[9] + Gx1[185]*Gu1[13] + Gx1[186]*Gu1[17] + Gx1[187]*Gu1[21] + Gx1[188]*Gu1[25] + Gx1[189]*Gu1[29] + Gx1[190]*Gu1[33] + Gx1[191]*Gu1[37] + Gx1[192]*Gu1[41] + Gx1[193]*Gu1[45] + Gx1[194]*Gu1[49] + Gx1[195]*Gu1[53];
Gu2[54] = + Gx1[182]*Gu1[2] + Gx1[183]*Gu1[6] + Gx1[184]*Gu1[10] + Gx1[185]*Gu1[14] + Gx1[186]*Gu1[18] + Gx1[187]*Gu1[22] + Gx1[188]*Gu1[26] + Gx1[189]*Gu1[30] + Gx1[190]*Gu1[34] + Gx1[191]*Gu1[38] + Gx1[192]*Gu1[42] + Gx1[193]*Gu1[46] + Gx1[194]*Gu1[50] + Gx1[195]*Gu1[54];
Gu2[55] = + Gx1[182]*Gu1[3] + Gx1[183]*Gu1[7] + Gx1[184]*Gu1[11] + Gx1[185]*Gu1[15] + Gx1[186]*Gu1[19] + Gx1[187]*Gu1[23] + Gx1[188]*Gu1[27] + Gx1[189]*Gu1[31] + Gx1[190]*Gu1[35] + Gx1[191]*Gu1[39] + Gx1[192]*Gu1[43] + Gx1[193]*Gu1[47] + Gx1[194]*Gu1[51] + Gx1[195]*Gu1[55];
}

void nmpc_moveGuE( real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = Gu1[0];
Gu2[1] = Gu1[1];
Gu2[2] = Gu1[2];
Gu2[3] = Gu1[3];
Gu2[4] = Gu1[4];
Gu2[5] = Gu1[5];
Gu2[6] = Gu1[6];
Gu2[7] = Gu1[7];
Gu2[8] = Gu1[8];
Gu2[9] = Gu1[9];
Gu2[10] = Gu1[10];
Gu2[11] = Gu1[11];
Gu2[12] = Gu1[12];
Gu2[13] = Gu1[13];
Gu2[14] = Gu1[14];
Gu2[15] = Gu1[15];
Gu2[16] = Gu1[16];
Gu2[17] = Gu1[17];
Gu2[18] = Gu1[18];
Gu2[19] = Gu1[19];
Gu2[20] = Gu1[20];
Gu2[21] = Gu1[21];
Gu2[22] = Gu1[22];
Gu2[23] = Gu1[23];
Gu2[24] = Gu1[24];
Gu2[25] = Gu1[25];
Gu2[26] = Gu1[26];
Gu2[27] = Gu1[27];
Gu2[28] = Gu1[28];
Gu2[29] = Gu1[29];
Gu2[30] = Gu1[30];
Gu2[31] = Gu1[31];
Gu2[32] = Gu1[32];
Gu2[33] = Gu1[33];
Gu2[34] = Gu1[34];
Gu2[35] = Gu1[35];
Gu2[36] = Gu1[36];
Gu2[37] = Gu1[37];
Gu2[38] = Gu1[38];
Gu2[39] = Gu1[39];
Gu2[40] = Gu1[40];
Gu2[41] = Gu1[41];
Gu2[42] = Gu1[42];
Gu2[43] = Gu1[43];
Gu2[44] = Gu1[44];
Gu2[45] = Gu1[45];
Gu2[46] = Gu1[46];
Gu2[47] = Gu1[47];
Gu2[48] = Gu1[48];
Gu2[49] = Gu1[49];
Gu2[50] = Gu1[50];
Gu2[51] = Gu1[51];
Gu2[52] = Gu1[52];
Gu2[53] = Gu1[53];
Gu2[54] = Gu1[54];
Gu2[55] = Gu1[55];
}

void nmpc_setBlockH11( int iRow, int iCol, real_t* const Gu1, real_t* const Gu2 )
{
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 14)] += + Gu1[0]*Gu2[0] + Gu1[4]*Gu2[4] + Gu1[8]*Gu2[8] + Gu1[12]*Gu2[12] + Gu1[16]*Gu2[16] + Gu1[20]*Gu2[20] + Gu1[24]*Gu2[24] + Gu1[28]*Gu2[28] + Gu1[32]*Gu2[32] + Gu1[36]*Gu2[36] + Gu1[40]*Gu2[40] + Gu1[44]*Gu2[44] + Gu1[48]*Gu2[48] + Gu1[52]*Gu2[52];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 15)] += + Gu1[0]*Gu2[1] + Gu1[4]*Gu2[5] + Gu1[8]*Gu2[9] + Gu1[12]*Gu2[13] + Gu1[16]*Gu2[17] + Gu1[20]*Gu2[21] + Gu1[24]*Gu2[25] + Gu1[28]*Gu2[29] + Gu1[32]*Gu2[33] + Gu1[36]*Gu2[37] + Gu1[40]*Gu2[41] + Gu1[44]*Gu2[45] + Gu1[48]*Gu2[49] + Gu1[52]*Gu2[53];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 16)] += + Gu1[0]*Gu2[2] + Gu1[4]*Gu2[6] + Gu1[8]*Gu2[10] + Gu1[12]*Gu2[14] + Gu1[16]*Gu2[18] + Gu1[20]*Gu2[22] + Gu1[24]*Gu2[26] + Gu1[28]*Gu2[30] + Gu1[32]*Gu2[34] + Gu1[36]*Gu2[38] + Gu1[40]*Gu2[42] + Gu1[44]*Gu2[46] + Gu1[48]*Gu2[50] + Gu1[52]*Gu2[54];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 17)] += + Gu1[0]*Gu2[3] + Gu1[4]*Gu2[7] + Gu1[8]*Gu2[11] + Gu1[12]*Gu2[15] + Gu1[16]*Gu2[19] + Gu1[20]*Gu2[23] + Gu1[24]*Gu2[27] + Gu1[28]*Gu2[31] + Gu1[32]*Gu2[35] + Gu1[36]*Gu2[39] + Gu1[40]*Gu2[43] + Gu1[44]*Gu2[47] + Gu1[48]*Gu2[51] + Gu1[52]*Gu2[55];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 14)] += + Gu1[1]*Gu2[0] + Gu1[5]*Gu2[4] + Gu1[9]*Gu2[8] + Gu1[13]*Gu2[12] + Gu1[17]*Gu2[16] + Gu1[21]*Gu2[20] + Gu1[25]*Gu2[24] + Gu1[29]*Gu2[28] + Gu1[33]*Gu2[32] + Gu1[37]*Gu2[36] + Gu1[41]*Gu2[40] + Gu1[45]*Gu2[44] + Gu1[49]*Gu2[48] + Gu1[53]*Gu2[52];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 15)] += + Gu1[1]*Gu2[1] + Gu1[5]*Gu2[5] + Gu1[9]*Gu2[9] + Gu1[13]*Gu2[13] + Gu1[17]*Gu2[17] + Gu1[21]*Gu2[21] + Gu1[25]*Gu2[25] + Gu1[29]*Gu2[29] + Gu1[33]*Gu2[33] + Gu1[37]*Gu2[37] + Gu1[41]*Gu2[41] + Gu1[45]*Gu2[45] + Gu1[49]*Gu2[49] + Gu1[53]*Gu2[53];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 16)] += + Gu1[1]*Gu2[2] + Gu1[5]*Gu2[6] + Gu1[9]*Gu2[10] + Gu1[13]*Gu2[14] + Gu1[17]*Gu2[18] + Gu1[21]*Gu2[22] + Gu1[25]*Gu2[26] + Gu1[29]*Gu2[30] + Gu1[33]*Gu2[34] + Gu1[37]*Gu2[38] + Gu1[41]*Gu2[42] + Gu1[45]*Gu2[46] + Gu1[49]*Gu2[50] + Gu1[53]*Gu2[54];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 17)] += + Gu1[1]*Gu2[3] + Gu1[5]*Gu2[7] + Gu1[9]*Gu2[11] + Gu1[13]*Gu2[15] + Gu1[17]*Gu2[19] + Gu1[21]*Gu2[23] + Gu1[25]*Gu2[27] + Gu1[29]*Gu2[31] + Gu1[33]*Gu2[35] + Gu1[37]*Gu2[39] + Gu1[41]*Gu2[43] + Gu1[45]*Gu2[47] + Gu1[49]*Gu2[51] + Gu1[53]*Gu2[55];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 14)] += + Gu1[2]*Gu2[0] + Gu1[6]*Gu2[4] + Gu1[10]*Gu2[8] + Gu1[14]*Gu2[12] + Gu1[18]*Gu2[16] + Gu1[22]*Gu2[20] + Gu1[26]*Gu2[24] + Gu1[30]*Gu2[28] + Gu1[34]*Gu2[32] + Gu1[38]*Gu2[36] + Gu1[42]*Gu2[40] + Gu1[46]*Gu2[44] + Gu1[50]*Gu2[48] + Gu1[54]*Gu2[52];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 15)] += + Gu1[2]*Gu2[1] + Gu1[6]*Gu2[5] + Gu1[10]*Gu2[9] + Gu1[14]*Gu2[13] + Gu1[18]*Gu2[17] + Gu1[22]*Gu2[21] + Gu1[26]*Gu2[25] + Gu1[30]*Gu2[29] + Gu1[34]*Gu2[33] + Gu1[38]*Gu2[37] + Gu1[42]*Gu2[41] + Gu1[46]*Gu2[45] + Gu1[50]*Gu2[49] + Gu1[54]*Gu2[53];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 16)] += + Gu1[2]*Gu2[2] + Gu1[6]*Gu2[6] + Gu1[10]*Gu2[10] + Gu1[14]*Gu2[14] + Gu1[18]*Gu2[18] + Gu1[22]*Gu2[22] + Gu1[26]*Gu2[26] + Gu1[30]*Gu2[30] + Gu1[34]*Gu2[34] + Gu1[38]*Gu2[38] + Gu1[42]*Gu2[42] + Gu1[46]*Gu2[46] + Gu1[50]*Gu2[50] + Gu1[54]*Gu2[54];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 17)] += + Gu1[2]*Gu2[3] + Gu1[6]*Gu2[7] + Gu1[10]*Gu2[11] + Gu1[14]*Gu2[15] + Gu1[18]*Gu2[19] + Gu1[22]*Gu2[23] + Gu1[26]*Gu2[27] + Gu1[30]*Gu2[31] + Gu1[34]*Gu2[35] + Gu1[38]*Gu2[39] + Gu1[42]*Gu2[43] + Gu1[46]*Gu2[47] + Gu1[50]*Gu2[51] + Gu1[54]*Gu2[55];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 14)] += + Gu1[3]*Gu2[0] + Gu1[7]*Gu2[4] + Gu1[11]*Gu2[8] + Gu1[15]*Gu2[12] + Gu1[19]*Gu2[16] + Gu1[23]*Gu2[20] + Gu1[27]*Gu2[24] + Gu1[31]*Gu2[28] + Gu1[35]*Gu2[32] + Gu1[39]*Gu2[36] + Gu1[43]*Gu2[40] + Gu1[47]*Gu2[44] + Gu1[51]*Gu2[48] + Gu1[55]*Gu2[52];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 15)] += + Gu1[3]*Gu2[1] + Gu1[7]*Gu2[5] + Gu1[11]*Gu2[9] + Gu1[15]*Gu2[13] + Gu1[19]*Gu2[17] + Gu1[23]*Gu2[21] + Gu1[27]*Gu2[25] + Gu1[31]*Gu2[29] + Gu1[35]*Gu2[33] + Gu1[39]*Gu2[37] + Gu1[43]*Gu2[41] + Gu1[47]*Gu2[45] + Gu1[51]*Gu2[49] + Gu1[55]*Gu2[53];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 16)] += + Gu1[3]*Gu2[2] + Gu1[7]*Gu2[6] + Gu1[11]*Gu2[10] + Gu1[15]*Gu2[14] + Gu1[19]*Gu2[18] + Gu1[23]*Gu2[22] + Gu1[27]*Gu2[26] + Gu1[31]*Gu2[30] + Gu1[35]*Gu2[34] + Gu1[39]*Gu2[38] + Gu1[43]*Gu2[42] + Gu1[47]*Gu2[46] + Gu1[51]*Gu2[50] + Gu1[55]*Gu2[54];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 17)] += + Gu1[3]*Gu2[3] + Gu1[7]*Gu2[7] + Gu1[11]*Gu2[11] + Gu1[15]*Gu2[15] + Gu1[19]*Gu2[19] + Gu1[23]*Gu2[23] + Gu1[27]*Gu2[27] + Gu1[31]*Gu2[31] + Gu1[35]*Gu2[35] + Gu1[39]*Gu2[39] + Gu1[43]*Gu2[43] + Gu1[47]*Gu2[47] + Gu1[51]*Gu2[51] + Gu1[55]*Gu2[55];
}

void nmpc_setBlockH11_R1( int iRow, int iCol, real_t* const R11 )
{
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 14)] = R11[0];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 15)] = R11[1];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 16)] = R11[2];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 17)] = R11[3];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 14)] = R11[4];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 15)] = R11[5];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 16)] = R11[6];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 17)] = R11[7];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 14)] = R11[8];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 15)] = R11[9];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 16)] = R11[10];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 17)] = R11[11];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 14)] = R11[12];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 15)] = R11[13];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 16)] = R11[14];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 17)] = R11[15];
}

void nmpc_zeroBlockH11( int iRow, int iCol )
{
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 14)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 15)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 16)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 17)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 14)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 15)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 16)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 17)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 14)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 15)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 16)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 17)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 14)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 15)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 16)] = 0.0000000000000000e+00;
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 17)] = 0.0000000000000000e+00;
}

void nmpc_copyHTH( int iRow, int iCol )
{
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 14)] = nmpcWorkspace.H[(iCol * 856 + 2996) + (iRow * 4 + 14)];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 15)] = nmpcWorkspace.H[(iCol * 856 + 3210) + (iRow * 4 + 14)];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 16)] = nmpcWorkspace.H[(iCol * 856 + 3424) + (iRow * 4 + 14)];
nmpcWorkspace.H[(iRow * 856 + 2996) + (iCol * 4 + 17)] = nmpcWorkspace.H[(iCol * 856 + 3638) + (iRow * 4 + 14)];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 14)] = nmpcWorkspace.H[(iCol * 856 + 2996) + (iRow * 4 + 15)];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 15)] = nmpcWorkspace.H[(iCol * 856 + 3210) + (iRow * 4 + 15)];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 16)] = nmpcWorkspace.H[(iCol * 856 + 3424) + (iRow * 4 + 15)];
nmpcWorkspace.H[(iRow * 856 + 3210) + (iCol * 4 + 17)] = nmpcWorkspace.H[(iCol * 856 + 3638) + (iRow * 4 + 15)];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 14)] = nmpcWorkspace.H[(iCol * 856 + 2996) + (iRow * 4 + 16)];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 15)] = nmpcWorkspace.H[(iCol * 856 + 3210) + (iRow * 4 + 16)];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 16)] = nmpcWorkspace.H[(iCol * 856 + 3424) + (iRow * 4 + 16)];
nmpcWorkspace.H[(iRow * 856 + 3424) + (iCol * 4 + 17)] = nmpcWorkspace.H[(iCol * 856 + 3638) + (iRow * 4 + 16)];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 14)] = nmpcWorkspace.H[(iCol * 856 + 2996) + (iRow * 4 + 17)];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 15)] = nmpcWorkspace.H[(iCol * 856 + 3210) + (iRow * 4 + 17)];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 16)] = nmpcWorkspace.H[(iCol * 856 + 3424) + (iRow * 4 + 17)];
nmpcWorkspace.H[(iRow * 856 + 3638) + (iCol * 4 + 17)] = nmpcWorkspace.H[(iCol * 856 + 3638) + (iRow * 4 + 17)];
}

void nmpc_multQ1d( real_t* const Gx1, real_t* const dOld, real_t* const dNew )
{
dNew[0] = + Gx1[0]*dOld[0] + Gx1[1]*dOld[1] + Gx1[2]*dOld[2] + Gx1[3]*dOld[3] + Gx1[4]*dOld[4] + Gx1[5]*dOld[5] + Gx1[6]*dOld[6] + Gx1[7]*dOld[7] + Gx1[8]*dOld[8] + Gx1[9]*dOld[9] + Gx1[10]*dOld[10] + Gx1[11]*dOld[11] + Gx1[12]*dOld[12] + Gx1[13]*dOld[13];
dNew[1] = + Gx1[14]*dOld[0] + Gx1[15]*dOld[1] + Gx1[16]*dOld[2] + Gx1[17]*dOld[3] + Gx1[18]*dOld[4] + Gx1[19]*dOld[5] + Gx1[20]*dOld[6] + Gx1[21]*dOld[7] + Gx1[22]*dOld[8] + Gx1[23]*dOld[9] + Gx1[24]*dOld[10] + Gx1[25]*dOld[11] + Gx1[26]*dOld[12] + Gx1[27]*dOld[13];
dNew[2] = + Gx1[28]*dOld[0] + Gx1[29]*dOld[1] + Gx1[30]*dOld[2] + Gx1[31]*dOld[3] + Gx1[32]*dOld[4] + Gx1[33]*dOld[5] + Gx1[34]*dOld[6] + Gx1[35]*dOld[7] + Gx1[36]*dOld[8] + Gx1[37]*dOld[9] + Gx1[38]*dOld[10] + Gx1[39]*dOld[11] + Gx1[40]*dOld[12] + Gx1[41]*dOld[13];
dNew[3] = + Gx1[42]*dOld[0] + Gx1[43]*dOld[1] + Gx1[44]*dOld[2] + Gx1[45]*dOld[3] + Gx1[46]*dOld[4] + Gx1[47]*dOld[5] + Gx1[48]*dOld[6] + Gx1[49]*dOld[7] + Gx1[50]*dOld[8] + Gx1[51]*dOld[9] + Gx1[52]*dOld[10] + Gx1[53]*dOld[11] + Gx1[54]*dOld[12] + Gx1[55]*dOld[13];
dNew[4] = + Gx1[56]*dOld[0] + Gx1[57]*dOld[1] + Gx1[58]*dOld[2] + Gx1[59]*dOld[3] + Gx1[60]*dOld[4] + Gx1[61]*dOld[5] + Gx1[62]*dOld[6] + Gx1[63]*dOld[7] + Gx1[64]*dOld[8] + Gx1[65]*dOld[9] + Gx1[66]*dOld[10] + Gx1[67]*dOld[11] + Gx1[68]*dOld[12] + Gx1[69]*dOld[13];
dNew[5] = + Gx1[70]*dOld[0] + Gx1[71]*dOld[1] + Gx1[72]*dOld[2] + Gx1[73]*dOld[3] + Gx1[74]*dOld[4] + Gx1[75]*dOld[5] + Gx1[76]*dOld[6] + Gx1[77]*dOld[7] + Gx1[78]*dOld[8] + Gx1[79]*dOld[9] + Gx1[80]*dOld[10] + Gx1[81]*dOld[11] + Gx1[82]*dOld[12] + Gx1[83]*dOld[13];
dNew[6] = + Gx1[84]*dOld[0] + Gx1[85]*dOld[1] + Gx1[86]*dOld[2] + Gx1[87]*dOld[3] + Gx1[88]*dOld[4] + Gx1[89]*dOld[5] + Gx1[90]*dOld[6] + Gx1[91]*dOld[7] + Gx1[92]*dOld[8] + Gx1[93]*dOld[9] + Gx1[94]*dOld[10] + Gx1[95]*dOld[11] + Gx1[96]*dOld[12] + Gx1[97]*dOld[13];
dNew[7] = + Gx1[98]*dOld[0] + Gx1[99]*dOld[1] + Gx1[100]*dOld[2] + Gx1[101]*dOld[3] + Gx1[102]*dOld[4] + Gx1[103]*dOld[5] + Gx1[104]*dOld[6] + Gx1[105]*dOld[7] + Gx1[106]*dOld[8] + Gx1[107]*dOld[9] + Gx1[108]*dOld[10] + Gx1[109]*dOld[11] + Gx1[110]*dOld[12] + Gx1[111]*dOld[13];
dNew[8] = + Gx1[112]*dOld[0] + Gx1[113]*dOld[1] + Gx1[114]*dOld[2] + Gx1[115]*dOld[3] + Gx1[116]*dOld[4] + Gx1[117]*dOld[5] + Gx1[118]*dOld[6] + Gx1[119]*dOld[7] + Gx1[120]*dOld[8] + Gx1[121]*dOld[9] + Gx1[122]*dOld[10] + Gx1[123]*dOld[11] + Gx1[124]*dOld[12] + Gx1[125]*dOld[13];
dNew[9] = + Gx1[126]*dOld[0] + Gx1[127]*dOld[1] + Gx1[128]*dOld[2] + Gx1[129]*dOld[3] + Gx1[130]*dOld[4] + Gx1[131]*dOld[5] + Gx1[132]*dOld[6] + Gx1[133]*dOld[7] + Gx1[134]*dOld[8] + Gx1[135]*dOld[9] + Gx1[136]*dOld[10] + Gx1[137]*dOld[11] + Gx1[138]*dOld[12] + Gx1[139]*dOld[13];
dNew[10] = + Gx1[140]*dOld[0] + Gx1[141]*dOld[1] + Gx1[142]*dOld[2] + Gx1[143]*dOld[3] + Gx1[144]*dOld[4] + Gx1[145]*dOld[5] + Gx1[146]*dOld[6] + Gx1[147]*dOld[7] + Gx1[148]*dOld[8] + Gx1[149]*dOld[9] + Gx1[150]*dOld[10] + Gx1[151]*dOld[11] + Gx1[152]*dOld[12] + Gx1[153]*dOld[13];
dNew[11] = + Gx1[154]*dOld[0] + Gx1[155]*dOld[1] + Gx1[156]*dOld[2] + Gx1[157]*dOld[3] + Gx1[158]*dOld[4] + Gx1[159]*dOld[5] + Gx1[160]*dOld[6] + Gx1[161]*dOld[7] + Gx1[162]*dOld[8] + Gx1[163]*dOld[9] + Gx1[164]*dOld[10] + Gx1[165]*dOld[11] + Gx1[166]*dOld[12] + Gx1[167]*dOld[13];
dNew[12] = + Gx1[168]*dOld[0] + Gx1[169]*dOld[1] + Gx1[170]*dOld[2] + Gx1[171]*dOld[3] + Gx1[172]*dOld[4] + Gx1[173]*dOld[5] + Gx1[174]*dOld[6] + Gx1[175]*dOld[7] + Gx1[176]*dOld[8] + Gx1[177]*dOld[9] + Gx1[178]*dOld[10] + Gx1[179]*dOld[11] + Gx1[180]*dOld[12] + Gx1[181]*dOld[13];
dNew[13] = + Gx1[182]*dOld[0] + Gx1[183]*dOld[1] + Gx1[184]*dOld[2] + Gx1[185]*dOld[3] + Gx1[186]*dOld[4] + Gx1[187]*dOld[5] + Gx1[188]*dOld[6] + Gx1[189]*dOld[7] + Gx1[190]*dOld[8] + Gx1[191]*dOld[9] + Gx1[192]*dOld[10] + Gx1[193]*dOld[11] + Gx1[194]*dOld[12] + Gx1[195]*dOld[13];
}

void nmpc_multQN1d( real_t* const QN1, real_t* const dOld, real_t* const dNew )
{
dNew[0] = + nmpcWorkspace.QN1[0]*dOld[0] + nmpcWorkspace.QN1[1]*dOld[1] + nmpcWorkspace.QN1[2]*dOld[2] + nmpcWorkspace.QN1[3]*dOld[3] + nmpcWorkspace.QN1[4]*dOld[4] + nmpcWorkspace.QN1[5]*dOld[5] + nmpcWorkspace.QN1[6]*dOld[6] + nmpcWorkspace.QN1[7]*dOld[7] + nmpcWorkspace.QN1[8]*dOld[8] + nmpcWorkspace.QN1[9]*dOld[9] + nmpcWorkspace.QN1[10]*dOld[10] + nmpcWorkspace.QN1[11]*dOld[11] + nmpcWorkspace.QN1[12]*dOld[12] + nmpcWorkspace.QN1[13]*dOld[13];
dNew[1] = + nmpcWorkspace.QN1[14]*dOld[0] + nmpcWorkspace.QN1[15]*dOld[1] + nmpcWorkspace.QN1[16]*dOld[2] + nmpcWorkspace.QN1[17]*dOld[3] + nmpcWorkspace.QN1[18]*dOld[4] + nmpcWorkspace.QN1[19]*dOld[5] + nmpcWorkspace.QN1[20]*dOld[6] + nmpcWorkspace.QN1[21]*dOld[7] + nmpcWorkspace.QN1[22]*dOld[8] + nmpcWorkspace.QN1[23]*dOld[9] + nmpcWorkspace.QN1[24]*dOld[10] + nmpcWorkspace.QN1[25]*dOld[11] + nmpcWorkspace.QN1[26]*dOld[12] + nmpcWorkspace.QN1[27]*dOld[13];
dNew[2] = + nmpcWorkspace.QN1[28]*dOld[0] + nmpcWorkspace.QN1[29]*dOld[1] + nmpcWorkspace.QN1[30]*dOld[2] + nmpcWorkspace.QN1[31]*dOld[3] + nmpcWorkspace.QN1[32]*dOld[4] + nmpcWorkspace.QN1[33]*dOld[5] + nmpcWorkspace.QN1[34]*dOld[6] + nmpcWorkspace.QN1[35]*dOld[7] + nmpcWorkspace.QN1[36]*dOld[8] + nmpcWorkspace.QN1[37]*dOld[9] + nmpcWorkspace.QN1[38]*dOld[10] + nmpcWorkspace.QN1[39]*dOld[11] + nmpcWorkspace.QN1[40]*dOld[12] + nmpcWorkspace.QN1[41]*dOld[13];
dNew[3] = + nmpcWorkspace.QN1[42]*dOld[0] + nmpcWorkspace.QN1[43]*dOld[1] + nmpcWorkspace.QN1[44]*dOld[2] + nmpcWorkspace.QN1[45]*dOld[3] + nmpcWorkspace.QN1[46]*dOld[4] + nmpcWorkspace.QN1[47]*dOld[5] + nmpcWorkspace.QN1[48]*dOld[6] + nmpcWorkspace.QN1[49]*dOld[7] + nmpcWorkspace.QN1[50]*dOld[8] + nmpcWorkspace.QN1[51]*dOld[9] + nmpcWorkspace.QN1[52]*dOld[10] + nmpcWorkspace.QN1[53]*dOld[11] + nmpcWorkspace.QN1[54]*dOld[12] + nmpcWorkspace.QN1[55]*dOld[13];
dNew[4] = + nmpcWorkspace.QN1[56]*dOld[0] + nmpcWorkspace.QN1[57]*dOld[1] + nmpcWorkspace.QN1[58]*dOld[2] + nmpcWorkspace.QN1[59]*dOld[3] + nmpcWorkspace.QN1[60]*dOld[4] + nmpcWorkspace.QN1[61]*dOld[5] + nmpcWorkspace.QN1[62]*dOld[6] + nmpcWorkspace.QN1[63]*dOld[7] + nmpcWorkspace.QN1[64]*dOld[8] + nmpcWorkspace.QN1[65]*dOld[9] + nmpcWorkspace.QN1[66]*dOld[10] + nmpcWorkspace.QN1[67]*dOld[11] + nmpcWorkspace.QN1[68]*dOld[12] + nmpcWorkspace.QN1[69]*dOld[13];
dNew[5] = + nmpcWorkspace.QN1[70]*dOld[0] + nmpcWorkspace.QN1[71]*dOld[1] + nmpcWorkspace.QN1[72]*dOld[2] + nmpcWorkspace.QN1[73]*dOld[3] + nmpcWorkspace.QN1[74]*dOld[4] + nmpcWorkspace.QN1[75]*dOld[5] + nmpcWorkspace.QN1[76]*dOld[6] + nmpcWorkspace.QN1[77]*dOld[7] + nmpcWorkspace.QN1[78]*dOld[8] + nmpcWorkspace.QN1[79]*dOld[9] + nmpcWorkspace.QN1[80]*dOld[10] + nmpcWorkspace.QN1[81]*dOld[11] + nmpcWorkspace.QN1[82]*dOld[12] + nmpcWorkspace.QN1[83]*dOld[13];
dNew[6] = + nmpcWorkspace.QN1[84]*dOld[0] + nmpcWorkspace.QN1[85]*dOld[1] + nmpcWorkspace.QN1[86]*dOld[2] + nmpcWorkspace.QN1[87]*dOld[3] + nmpcWorkspace.QN1[88]*dOld[4] + nmpcWorkspace.QN1[89]*dOld[5] + nmpcWorkspace.QN1[90]*dOld[6] + nmpcWorkspace.QN1[91]*dOld[7] + nmpcWorkspace.QN1[92]*dOld[8] + nmpcWorkspace.QN1[93]*dOld[9] + nmpcWorkspace.QN1[94]*dOld[10] + nmpcWorkspace.QN1[95]*dOld[11] + nmpcWorkspace.QN1[96]*dOld[12] + nmpcWorkspace.QN1[97]*dOld[13];
dNew[7] = + nmpcWorkspace.QN1[98]*dOld[0] + nmpcWorkspace.QN1[99]*dOld[1] + nmpcWorkspace.QN1[100]*dOld[2] + nmpcWorkspace.QN1[101]*dOld[3] + nmpcWorkspace.QN1[102]*dOld[4] + nmpcWorkspace.QN1[103]*dOld[5] + nmpcWorkspace.QN1[104]*dOld[6] + nmpcWorkspace.QN1[105]*dOld[7] + nmpcWorkspace.QN1[106]*dOld[8] + nmpcWorkspace.QN1[107]*dOld[9] + nmpcWorkspace.QN1[108]*dOld[10] + nmpcWorkspace.QN1[109]*dOld[11] + nmpcWorkspace.QN1[110]*dOld[12] + nmpcWorkspace.QN1[111]*dOld[13];
dNew[8] = + nmpcWorkspace.QN1[112]*dOld[0] + nmpcWorkspace.QN1[113]*dOld[1] + nmpcWorkspace.QN1[114]*dOld[2] + nmpcWorkspace.QN1[115]*dOld[3] + nmpcWorkspace.QN1[116]*dOld[4] + nmpcWorkspace.QN1[117]*dOld[5] + nmpcWorkspace.QN1[118]*dOld[6] + nmpcWorkspace.QN1[119]*dOld[7] + nmpcWorkspace.QN1[120]*dOld[8] + nmpcWorkspace.QN1[121]*dOld[9] + nmpcWorkspace.QN1[122]*dOld[10] + nmpcWorkspace.QN1[123]*dOld[11] + nmpcWorkspace.QN1[124]*dOld[12] + nmpcWorkspace.QN1[125]*dOld[13];
dNew[9] = + nmpcWorkspace.QN1[126]*dOld[0] + nmpcWorkspace.QN1[127]*dOld[1] + nmpcWorkspace.QN1[128]*dOld[2] + nmpcWorkspace.QN1[129]*dOld[3] + nmpcWorkspace.QN1[130]*dOld[4] + nmpcWorkspace.QN1[131]*dOld[5] + nmpcWorkspace.QN1[132]*dOld[6] + nmpcWorkspace.QN1[133]*dOld[7] + nmpcWorkspace.QN1[134]*dOld[8] + nmpcWorkspace.QN1[135]*dOld[9] + nmpcWorkspace.QN1[136]*dOld[10] + nmpcWorkspace.QN1[137]*dOld[11] + nmpcWorkspace.QN1[138]*dOld[12] + nmpcWorkspace.QN1[139]*dOld[13];
dNew[10] = + nmpcWorkspace.QN1[140]*dOld[0] + nmpcWorkspace.QN1[141]*dOld[1] + nmpcWorkspace.QN1[142]*dOld[2] + nmpcWorkspace.QN1[143]*dOld[3] + nmpcWorkspace.QN1[144]*dOld[4] + nmpcWorkspace.QN1[145]*dOld[5] + nmpcWorkspace.QN1[146]*dOld[6] + nmpcWorkspace.QN1[147]*dOld[7] + nmpcWorkspace.QN1[148]*dOld[8] + nmpcWorkspace.QN1[149]*dOld[9] + nmpcWorkspace.QN1[150]*dOld[10] + nmpcWorkspace.QN1[151]*dOld[11] + nmpcWorkspace.QN1[152]*dOld[12] + nmpcWorkspace.QN1[153]*dOld[13];
dNew[11] = + nmpcWorkspace.QN1[154]*dOld[0] + nmpcWorkspace.QN1[155]*dOld[1] + nmpcWorkspace.QN1[156]*dOld[2] + nmpcWorkspace.QN1[157]*dOld[3] + nmpcWorkspace.QN1[158]*dOld[4] + nmpcWorkspace.QN1[159]*dOld[5] + nmpcWorkspace.QN1[160]*dOld[6] + nmpcWorkspace.QN1[161]*dOld[7] + nmpcWorkspace.QN1[162]*dOld[8] + nmpcWorkspace.QN1[163]*dOld[9] + nmpcWorkspace.QN1[164]*dOld[10] + nmpcWorkspace.QN1[165]*dOld[11] + nmpcWorkspace.QN1[166]*dOld[12] + nmpcWorkspace.QN1[167]*dOld[13];
dNew[12] = + nmpcWorkspace.QN1[168]*dOld[0] + nmpcWorkspace.QN1[169]*dOld[1] + nmpcWorkspace.QN1[170]*dOld[2] + nmpcWorkspace.QN1[171]*dOld[3] + nmpcWorkspace.QN1[172]*dOld[4] + nmpcWorkspace.QN1[173]*dOld[5] + nmpcWorkspace.QN1[174]*dOld[6] + nmpcWorkspace.QN1[175]*dOld[7] + nmpcWorkspace.QN1[176]*dOld[8] + nmpcWorkspace.QN1[177]*dOld[9] + nmpcWorkspace.QN1[178]*dOld[10] + nmpcWorkspace.QN1[179]*dOld[11] + nmpcWorkspace.QN1[180]*dOld[12] + nmpcWorkspace.QN1[181]*dOld[13];
dNew[13] = + nmpcWorkspace.QN1[182]*dOld[0] + nmpcWorkspace.QN1[183]*dOld[1] + nmpcWorkspace.QN1[184]*dOld[2] + nmpcWorkspace.QN1[185]*dOld[3] + nmpcWorkspace.QN1[186]*dOld[4] + nmpcWorkspace.QN1[187]*dOld[5] + nmpcWorkspace.QN1[188]*dOld[6] + nmpcWorkspace.QN1[189]*dOld[7] + nmpcWorkspace.QN1[190]*dOld[8] + nmpcWorkspace.QN1[191]*dOld[9] + nmpcWorkspace.QN1[192]*dOld[10] + nmpcWorkspace.QN1[193]*dOld[11] + nmpcWorkspace.QN1[194]*dOld[12] + nmpcWorkspace.QN1[195]*dOld[13];
}

void nmpc_multRDy( real_t* const R2, real_t* const Dy1, real_t* const RDy1 )
{
RDy1[0] = + R2[0]*Dy1[0] + R2[1]*Dy1[1] + R2[2]*Dy1[2] + R2[3]*Dy1[3] + R2[4]*Dy1[4] + R2[5]*Dy1[5] + R2[6]*Dy1[6] + R2[7]*Dy1[7] + R2[8]*Dy1[8] + R2[9]*Dy1[9] + R2[10]*Dy1[10] + R2[11]*Dy1[11] + R2[12]*Dy1[12];
RDy1[1] = + R2[13]*Dy1[0] + R2[14]*Dy1[1] + R2[15]*Dy1[2] + R2[16]*Dy1[3] + R2[17]*Dy1[4] + R2[18]*Dy1[5] + R2[19]*Dy1[6] + R2[20]*Dy1[7] + R2[21]*Dy1[8] + R2[22]*Dy1[9] + R2[23]*Dy1[10] + R2[24]*Dy1[11] + R2[25]*Dy1[12];
RDy1[2] = + R2[26]*Dy1[0] + R2[27]*Dy1[1] + R2[28]*Dy1[2] + R2[29]*Dy1[3] + R2[30]*Dy1[4] + R2[31]*Dy1[5] + R2[32]*Dy1[6] + R2[33]*Dy1[7] + R2[34]*Dy1[8] + R2[35]*Dy1[9] + R2[36]*Dy1[10] + R2[37]*Dy1[11] + R2[38]*Dy1[12];
RDy1[3] = + R2[39]*Dy1[0] + R2[40]*Dy1[1] + R2[41]*Dy1[2] + R2[42]*Dy1[3] + R2[43]*Dy1[4] + R2[44]*Dy1[5] + R2[45]*Dy1[6] + R2[46]*Dy1[7] + R2[47]*Dy1[8] + R2[48]*Dy1[9] + R2[49]*Dy1[10] + R2[50]*Dy1[11] + R2[51]*Dy1[12];
}

void nmpc_multQDy( real_t* const Q2, real_t* const Dy1, real_t* const QDy1 )
{
QDy1[0] = + Q2[0]*Dy1[0] + Q2[1]*Dy1[1] + Q2[2]*Dy1[2] + Q2[3]*Dy1[3] + Q2[4]*Dy1[4] + Q2[5]*Dy1[5] + Q2[6]*Dy1[6] + Q2[7]*Dy1[7] + Q2[8]*Dy1[8] + Q2[9]*Dy1[9] + Q2[10]*Dy1[10] + Q2[11]*Dy1[11] + Q2[12]*Dy1[12];
QDy1[1] = + Q2[13]*Dy1[0] + Q2[14]*Dy1[1] + Q2[15]*Dy1[2] + Q2[16]*Dy1[3] + Q2[17]*Dy1[4] + Q2[18]*Dy1[5] + Q2[19]*Dy1[6] + Q2[20]*Dy1[7] + Q2[21]*Dy1[8] + Q2[22]*Dy1[9] + Q2[23]*Dy1[10] + Q2[24]*Dy1[11] + Q2[25]*Dy1[12];
QDy1[2] = + Q2[26]*Dy1[0] + Q2[27]*Dy1[1] + Q2[28]*Dy1[2] + Q2[29]*Dy1[3] + Q2[30]*Dy1[4] + Q2[31]*Dy1[5] + Q2[32]*Dy1[6] + Q2[33]*Dy1[7] + Q2[34]*Dy1[8] + Q2[35]*Dy1[9] + Q2[36]*Dy1[10] + Q2[37]*Dy1[11] + Q2[38]*Dy1[12];
QDy1[3] = + Q2[39]*Dy1[0] + Q2[40]*Dy1[1] + Q2[41]*Dy1[2] + Q2[42]*Dy1[3] + Q2[43]*Dy1[4] + Q2[44]*Dy1[5] + Q2[45]*Dy1[6] + Q2[46]*Dy1[7] + Q2[47]*Dy1[8] + Q2[48]*Dy1[9] + Q2[49]*Dy1[10] + Q2[50]*Dy1[11] + Q2[51]*Dy1[12];
QDy1[4] = + Q2[52]*Dy1[0] + Q2[53]*Dy1[1] + Q2[54]*Dy1[2] + Q2[55]*Dy1[3] + Q2[56]*Dy1[4] + Q2[57]*Dy1[5] + Q2[58]*Dy1[6] + Q2[59]*Dy1[7] + Q2[60]*Dy1[8] + Q2[61]*Dy1[9] + Q2[62]*Dy1[10] + Q2[63]*Dy1[11] + Q2[64]*Dy1[12];
QDy1[5] = + Q2[65]*Dy1[0] + Q2[66]*Dy1[1] + Q2[67]*Dy1[2] + Q2[68]*Dy1[3] + Q2[69]*Dy1[4] + Q2[70]*Dy1[5] + Q2[71]*Dy1[6] + Q2[72]*Dy1[7] + Q2[73]*Dy1[8] + Q2[74]*Dy1[9] + Q2[75]*Dy1[10] + Q2[76]*Dy1[11] + Q2[77]*Dy1[12];
QDy1[6] = + Q2[78]*Dy1[0] + Q2[79]*Dy1[1] + Q2[80]*Dy1[2] + Q2[81]*Dy1[3] + Q2[82]*Dy1[4] + Q2[83]*Dy1[5] + Q2[84]*Dy1[6] + Q2[85]*Dy1[7] + Q2[86]*Dy1[8] + Q2[87]*Dy1[9] + Q2[88]*Dy1[10] + Q2[89]*Dy1[11] + Q2[90]*Dy1[12];
QDy1[7] = + Q2[91]*Dy1[0] + Q2[92]*Dy1[1] + Q2[93]*Dy1[2] + Q2[94]*Dy1[3] + Q2[95]*Dy1[4] + Q2[96]*Dy1[5] + Q2[97]*Dy1[6] + Q2[98]*Dy1[7] + Q2[99]*Dy1[8] + Q2[100]*Dy1[9] + Q2[101]*Dy1[10] + Q2[102]*Dy1[11] + Q2[103]*Dy1[12];
QDy1[8] = + Q2[104]*Dy1[0] + Q2[105]*Dy1[1] + Q2[106]*Dy1[2] + Q2[107]*Dy1[3] + Q2[108]*Dy1[4] + Q2[109]*Dy1[5] + Q2[110]*Dy1[6] + Q2[111]*Dy1[7] + Q2[112]*Dy1[8] + Q2[113]*Dy1[9] + Q2[114]*Dy1[10] + Q2[115]*Dy1[11] + Q2[116]*Dy1[12];
QDy1[9] = + Q2[117]*Dy1[0] + Q2[118]*Dy1[1] + Q2[119]*Dy1[2] + Q2[120]*Dy1[3] + Q2[121]*Dy1[4] + Q2[122]*Dy1[5] + Q2[123]*Dy1[6] + Q2[124]*Dy1[7] + Q2[125]*Dy1[8] + Q2[126]*Dy1[9] + Q2[127]*Dy1[10] + Q2[128]*Dy1[11] + Q2[129]*Dy1[12];
QDy1[10] = + Q2[130]*Dy1[0] + Q2[131]*Dy1[1] + Q2[132]*Dy1[2] + Q2[133]*Dy1[3] + Q2[134]*Dy1[4] + Q2[135]*Dy1[5] + Q2[136]*Dy1[6] + Q2[137]*Dy1[7] + Q2[138]*Dy1[8] + Q2[139]*Dy1[9] + Q2[140]*Dy1[10] + Q2[141]*Dy1[11] + Q2[142]*Dy1[12];
QDy1[11] = + Q2[143]*Dy1[0] + Q2[144]*Dy1[1] + Q2[145]*Dy1[2] + Q2[146]*Dy1[3] + Q2[147]*Dy1[4] + Q2[148]*Dy1[5] + Q2[149]*Dy1[6] + Q2[150]*Dy1[7] + Q2[151]*Dy1[8] + Q2[152]*Dy1[9] + Q2[153]*Dy1[10] + Q2[154]*Dy1[11] + Q2[155]*Dy1[12];
QDy1[12] = + Q2[156]*Dy1[0] + Q2[157]*Dy1[1] + Q2[158]*Dy1[2] + Q2[159]*Dy1[3] + Q2[160]*Dy1[4] + Q2[161]*Dy1[5] + Q2[162]*Dy1[6] + Q2[163]*Dy1[7] + Q2[164]*Dy1[8] + Q2[165]*Dy1[9] + Q2[166]*Dy1[10] + Q2[167]*Dy1[11] + Q2[168]*Dy1[12];
QDy1[13] = + Q2[169]*Dy1[0] + Q2[170]*Dy1[1] + Q2[171]*Dy1[2] + Q2[172]*Dy1[3] + Q2[173]*Dy1[4] + Q2[174]*Dy1[5] + Q2[175]*Dy1[6] + Q2[176]*Dy1[7] + Q2[177]*Dy1[8] + Q2[178]*Dy1[9] + Q2[179]*Dy1[10] + Q2[180]*Dy1[11] + Q2[181]*Dy1[12];
}

void nmpc_multEQDy( real_t* const E1, real_t* const QDy1, real_t* const U1 )
{
U1[0] += + E1[0]*QDy1[0] + E1[4]*QDy1[1] + E1[8]*QDy1[2] + E1[12]*QDy1[3] + E1[16]*QDy1[4] + E1[20]*QDy1[5] + E1[24]*QDy1[6] + E1[28]*QDy1[7] + E1[32]*QDy1[8] + E1[36]*QDy1[9] + E1[40]*QDy1[10] + E1[44]*QDy1[11] + E1[48]*QDy1[12] + E1[52]*QDy1[13];
U1[1] += + E1[1]*QDy1[0] + E1[5]*QDy1[1] + E1[9]*QDy1[2] + E1[13]*QDy1[3] + E1[17]*QDy1[4] + E1[21]*QDy1[5] + E1[25]*QDy1[6] + E1[29]*QDy1[7] + E1[33]*QDy1[8] + E1[37]*QDy1[9] + E1[41]*QDy1[10] + E1[45]*QDy1[11] + E1[49]*QDy1[12] + E1[53]*QDy1[13];
U1[2] += + E1[2]*QDy1[0] + E1[6]*QDy1[1] + E1[10]*QDy1[2] + E1[14]*QDy1[3] + E1[18]*QDy1[4] + E1[22]*QDy1[5] + E1[26]*QDy1[6] + E1[30]*QDy1[7] + E1[34]*QDy1[8] + E1[38]*QDy1[9] + E1[42]*QDy1[10] + E1[46]*QDy1[11] + E1[50]*QDy1[12] + E1[54]*QDy1[13];
U1[3] += + E1[3]*QDy1[0] + E1[7]*QDy1[1] + E1[11]*QDy1[2] + E1[15]*QDy1[3] + E1[19]*QDy1[4] + E1[23]*QDy1[5] + E1[27]*QDy1[6] + E1[31]*QDy1[7] + E1[35]*QDy1[8] + E1[39]*QDy1[9] + E1[43]*QDy1[10] + E1[47]*QDy1[11] + E1[51]*QDy1[12] + E1[55]*QDy1[13];
}

void nmpc_multQETGx( real_t* const E1, real_t* const Gx1, real_t* const H101 )
{
H101[0] += + E1[0]*Gx1[0] + E1[4]*Gx1[14] + E1[8]*Gx1[28] + E1[12]*Gx1[42] + E1[16]*Gx1[56] + E1[20]*Gx1[70] + E1[24]*Gx1[84] + E1[28]*Gx1[98] + E1[32]*Gx1[112] + E1[36]*Gx1[126] + E1[40]*Gx1[140] + E1[44]*Gx1[154] + E1[48]*Gx1[168] + E1[52]*Gx1[182];
H101[1] += + E1[0]*Gx1[1] + E1[4]*Gx1[15] + E1[8]*Gx1[29] + E1[12]*Gx1[43] + E1[16]*Gx1[57] + E1[20]*Gx1[71] + E1[24]*Gx1[85] + E1[28]*Gx1[99] + E1[32]*Gx1[113] + E1[36]*Gx1[127] + E1[40]*Gx1[141] + E1[44]*Gx1[155] + E1[48]*Gx1[169] + E1[52]*Gx1[183];
H101[2] += + E1[0]*Gx1[2] + E1[4]*Gx1[16] + E1[8]*Gx1[30] + E1[12]*Gx1[44] + E1[16]*Gx1[58] + E1[20]*Gx1[72] + E1[24]*Gx1[86] + E1[28]*Gx1[100] + E1[32]*Gx1[114] + E1[36]*Gx1[128] + E1[40]*Gx1[142] + E1[44]*Gx1[156] + E1[48]*Gx1[170] + E1[52]*Gx1[184];
H101[3] += + E1[0]*Gx1[3] + E1[4]*Gx1[17] + E1[8]*Gx1[31] + E1[12]*Gx1[45] + E1[16]*Gx1[59] + E1[20]*Gx1[73] + E1[24]*Gx1[87] + E1[28]*Gx1[101] + E1[32]*Gx1[115] + E1[36]*Gx1[129] + E1[40]*Gx1[143] + E1[44]*Gx1[157] + E1[48]*Gx1[171] + E1[52]*Gx1[185];
H101[4] += + E1[0]*Gx1[4] + E1[4]*Gx1[18] + E1[8]*Gx1[32] + E1[12]*Gx1[46] + E1[16]*Gx1[60] + E1[20]*Gx1[74] + E1[24]*Gx1[88] + E1[28]*Gx1[102] + E1[32]*Gx1[116] + E1[36]*Gx1[130] + E1[40]*Gx1[144] + E1[44]*Gx1[158] + E1[48]*Gx1[172] + E1[52]*Gx1[186];
H101[5] += + E1[0]*Gx1[5] + E1[4]*Gx1[19] + E1[8]*Gx1[33] + E1[12]*Gx1[47] + E1[16]*Gx1[61] + E1[20]*Gx1[75] + E1[24]*Gx1[89] + E1[28]*Gx1[103] + E1[32]*Gx1[117] + E1[36]*Gx1[131] + E1[40]*Gx1[145] + E1[44]*Gx1[159] + E1[48]*Gx1[173] + E1[52]*Gx1[187];
H101[6] += + E1[0]*Gx1[6] + E1[4]*Gx1[20] + E1[8]*Gx1[34] + E1[12]*Gx1[48] + E1[16]*Gx1[62] + E1[20]*Gx1[76] + E1[24]*Gx1[90] + E1[28]*Gx1[104] + E1[32]*Gx1[118] + E1[36]*Gx1[132] + E1[40]*Gx1[146] + E1[44]*Gx1[160] + E1[48]*Gx1[174] + E1[52]*Gx1[188];
H101[7] += + E1[0]*Gx1[7] + E1[4]*Gx1[21] + E1[8]*Gx1[35] + E1[12]*Gx1[49] + E1[16]*Gx1[63] + E1[20]*Gx1[77] + E1[24]*Gx1[91] + E1[28]*Gx1[105] + E1[32]*Gx1[119] + E1[36]*Gx1[133] + E1[40]*Gx1[147] + E1[44]*Gx1[161] + E1[48]*Gx1[175] + E1[52]*Gx1[189];
H101[8] += + E1[0]*Gx1[8] + E1[4]*Gx1[22] + E1[8]*Gx1[36] + E1[12]*Gx1[50] + E1[16]*Gx1[64] + E1[20]*Gx1[78] + E1[24]*Gx1[92] + E1[28]*Gx1[106] + E1[32]*Gx1[120] + E1[36]*Gx1[134] + E1[40]*Gx1[148] + E1[44]*Gx1[162] + E1[48]*Gx1[176] + E1[52]*Gx1[190];
H101[9] += + E1[0]*Gx1[9] + E1[4]*Gx1[23] + E1[8]*Gx1[37] + E1[12]*Gx1[51] + E1[16]*Gx1[65] + E1[20]*Gx1[79] + E1[24]*Gx1[93] + E1[28]*Gx1[107] + E1[32]*Gx1[121] + E1[36]*Gx1[135] + E1[40]*Gx1[149] + E1[44]*Gx1[163] + E1[48]*Gx1[177] + E1[52]*Gx1[191];
H101[10] += + E1[0]*Gx1[10] + E1[4]*Gx1[24] + E1[8]*Gx1[38] + E1[12]*Gx1[52] + E1[16]*Gx1[66] + E1[20]*Gx1[80] + E1[24]*Gx1[94] + E1[28]*Gx1[108] + E1[32]*Gx1[122] + E1[36]*Gx1[136] + E1[40]*Gx1[150] + E1[44]*Gx1[164] + E1[48]*Gx1[178] + E1[52]*Gx1[192];
H101[11] += + E1[0]*Gx1[11] + E1[4]*Gx1[25] + E1[8]*Gx1[39] + E1[12]*Gx1[53] + E1[16]*Gx1[67] + E1[20]*Gx1[81] + E1[24]*Gx1[95] + E1[28]*Gx1[109] + E1[32]*Gx1[123] + E1[36]*Gx1[137] + E1[40]*Gx1[151] + E1[44]*Gx1[165] + E1[48]*Gx1[179] + E1[52]*Gx1[193];
H101[12] += + E1[0]*Gx1[12] + E1[4]*Gx1[26] + E1[8]*Gx1[40] + E1[12]*Gx1[54] + E1[16]*Gx1[68] + E1[20]*Gx1[82] + E1[24]*Gx1[96] + E1[28]*Gx1[110] + E1[32]*Gx1[124] + E1[36]*Gx1[138] + E1[40]*Gx1[152] + E1[44]*Gx1[166] + E1[48]*Gx1[180] + E1[52]*Gx1[194];
H101[13] += + E1[0]*Gx1[13] + E1[4]*Gx1[27] + E1[8]*Gx1[41] + E1[12]*Gx1[55] + E1[16]*Gx1[69] + E1[20]*Gx1[83] + E1[24]*Gx1[97] + E1[28]*Gx1[111] + E1[32]*Gx1[125] + E1[36]*Gx1[139] + E1[40]*Gx1[153] + E1[44]*Gx1[167] + E1[48]*Gx1[181] + E1[52]*Gx1[195];
H101[14] += + E1[1]*Gx1[0] + E1[5]*Gx1[14] + E1[9]*Gx1[28] + E1[13]*Gx1[42] + E1[17]*Gx1[56] + E1[21]*Gx1[70] + E1[25]*Gx1[84] + E1[29]*Gx1[98] + E1[33]*Gx1[112] + E1[37]*Gx1[126] + E1[41]*Gx1[140] + E1[45]*Gx1[154] + E1[49]*Gx1[168] + E1[53]*Gx1[182];
H101[15] += + E1[1]*Gx1[1] + E1[5]*Gx1[15] + E1[9]*Gx1[29] + E1[13]*Gx1[43] + E1[17]*Gx1[57] + E1[21]*Gx1[71] + E1[25]*Gx1[85] + E1[29]*Gx1[99] + E1[33]*Gx1[113] + E1[37]*Gx1[127] + E1[41]*Gx1[141] + E1[45]*Gx1[155] + E1[49]*Gx1[169] + E1[53]*Gx1[183];
H101[16] += + E1[1]*Gx1[2] + E1[5]*Gx1[16] + E1[9]*Gx1[30] + E1[13]*Gx1[44] + E1[17]*Gx1[58] + E1[21]*Gx1[72] + E1[25]*Gx1[86] + E1[29]*Gx1[100] + E1[33]*Gx1[114] + E1[37]*Gx1[128] + E1[41]*Gx1[142] + E1[45]*Gx1[156] + E1[49]*Gx1[170] + E1[53]*Gx1[184];
H101[17] += + E1[1]*Gx1[3] + E1[5]*Gx1[17] + E1[9]*Gx1[31] + E1[13]*Gx1[45] + E1[17]*Gx1[59] + E1[21]*Gx1[73] + E1[25]*Gx1[87] + E1[29]*Gx1[101] + E1[33]*Gx1[115] + E1[37]*Gx1[129] + E1[41]*Gx1[143] + E1[45]*Gx1[157] + E1[49]*Gx1[171] + E1[53]*Gx1[185];
H101[18] += + E1[1]*Gx1[4] + E1[5]*Gx1[18] + E1[9]*Gx1[32] + E1[13]*Gx1[46] + E1[17]*Gx1[60] + E1[21]*Gx1[74] + E1[25]*Gx1[88] + E1[29]*Gx1[102] + E1[33]*Gx1[116] + E1[37]*Gx1[130] + E1[41]*Gx1[144] + E1[45]*Gx1[158] + E1[49]*Gx1[172] + E1[53]*Gx1[186];
H101[19] += + E1[1]*Gx1[5] + E1[5]*Gx1[19] + E1[9]*Gx1[33] + E1[13]*Gx1[47] + E1[17]*Gx1[61] + E1[21]*Gx1[75] + E1[25]*Gx1[89] + E1[29]*Gx1[103] + E1[33]*Gx1[117] + E1[37]*Gx1[131] + E1[41]*Gx1[145] + E1[45]*Gx1[159] + E1[49]*Gx1[173] + E1[53]*Gx1[187];
H101[20] += + E1[1]*Gx1[6] + E1[5]*Gx1[20] + E1[9]*Gx1[34] + E1[13]*Gx1[48] + E1[17]*Gx1[62] + E1[21]*Gx1[76] + E1[25]*Gx1[90] + E1[29]*Gx1[104] + E1[33]*Gx1[118] + E1[37]*Gx1[132] + E1[41]*Gx1[146] + E1[45]*Gx1[160] + E1[49]*Gx1[174] + E1[53]*Gx1[188];
H101[21] += + E1[1]*Gx1[7] + E1[5]*Gx1[21] + E1[9]*Gx1[35] + E1[13]*Gx1[49] + E1[17]*Gx1[63] + E1[21]*Gx1[77] + E1[25]*Gx1[91] + E1[29]*Gx1[105] + E1[33]*Gx1[119] + E1[37]*Gx1[133] + E1[41]*Gx1[147] + E1[45]*Gx1[161] + E1[49]*Gx1[175] + E1[53]*Gx1[189];
H101[22] += + E1[1]*Gx1[8] + E1[5]*Gx1[22] + E1[9]*Gx1[36] + E1[13]*Gx1[50] + E1[17]*Gx1[64] + E1[21]*Gx1[78] + E1[25]*Gx1[92] + E1[29]*Gx1[106] + E1[33]*Gx1[120] + E1[37]*Gx1[134] + E1[41]*Gx1[148] + E1[45]*Gx1[162] + E1[49]*Gx1[176] + E1[53]*Gx1[190];
H101[23] += + E1[1]*Gx1[9] + E1[5]*Gx1[23] + E1[9]*Gx1[37] + E1[13]*Gx1[51] + E1[17]*Gx1[65] + E1[21]*Gx1[79] + E1[25]*Gx1[93] + E1[29]*Gx1[107] + E1[33]*Gx1[121] + E1[37]*Gx1[135] + E1[41]*Gx1[149] + E1[45]*Gx1[163] + E1[49]*Gx1[177] + E1[53]*Gx1[191];
H101[24] += + E1[1]*Gx1[10] + E1[5]*Gx1[24] + E1[9]*Gx1[38] + E1[13]*Gx1[52] + E1[17]*Gx1[66] + E1[21]*Gx1[80] + E1[25]*Gx1[94] + E1[29]*Gx1[108] + E1[33]*Gx1[122] + E1[37]*Gx1[136] + E1[41]*Gx1[150] + E1[45]*Gx1[164] + E1[49]*Gx1[178] + E1[53]*Gx1[192];
H101[25] += + E1[1]*Gx1[11] + E1[5]*Gx1[25] + E1[9]*Gx1[39] + E1[13]*Gx1[53] + E1[17]*Gx1[67] + E1[21]*Gx1[81] + E1[25]*Gx1[95] + E1[29]*Gx1[109] + E1[33]*Gx1[123] + E1[37]*Gx1[137] + E1[41]*Gx1[151] + E1[45]*Gx1[165] + E1[49]*Gx1[179] + E1[53]*Gx1[193];
H101[26] += + E1[1]*Gx1[12] + E1[5]*Gx1[26] + E1[9]*Gx1[40] + E1[13]*Gx1[54] + E1[17]*Gx1[68] + E1[21]*Gx1[82] + E1[25]*Gx1[96] + E1[29]*Gx1[110] + E1[33]*Gx1[124] + E1[37]*Gx1[138] + E1[41]*Gx1[152] + E1[45]*Gx1[166] + E1[49]*Gx1[180] + E1[53]*Gx1[194];
H101[27] += + E1[1]*Gx1[13] + E1[5]*Gx1[27] + E1[9]*Gx1[41] + E1[13]*Gx1[55] + E1[17]*Gx1[69] + E1[21]*Gx1[83] + E1[25]*Gx1[97] + E1[29]*Gx1[111] + E1[33]*Gx1[125] + E1[37]*Gx1[139] + E1[41]*Gx1[153] + E1[45]*Gx1[167] + E1[49]*Gx1[181] + E1[53]*Gx1[195];
H101[28] += + E1[2]*Gx1[0] + E1[6]*Gx1[14] + E1[10]*Gx1[28] + E1[14]*Gx1[42] + E1[18]*Gx1[56] + E1[22]*Gx1[70] + E1[26]*Gx1[84] + E1[30]*Gx1[98] + E1[34]*Gx1[112] + E1[38]*Gx1[126] + E1[42]*Gx1[140] + E1[46]*Gx1[154] + E1[50]*Gx1[168] + E1[54]*Gx1[182];
H101[29] += + E1[2]*Gx1[1] + E1[6]*Gx1[15] + E1[10]*Gx1[29] + E1[14]*Gx1[43] + E1[18]*Gx1[57] + E1[22]*Gx1[71] + E1[26]*Gx1[85] + E1[30]*Gx1[99] + E1[34]*Gx1[113] + E1[38]*Gx1[127] + E1[42]*Gx1[141] + E1[46]*Gx1[155] + E1[50]*Gx1[169] + E1[54]*Gx1[183];
H101[30] += + E1[2]*Gx1[2] + E1[6]*Gx1[16] + E1[10]*Gx1[30] + E1[14]*Gx1[44] + E1[18]*Gx1[58] + E1[22]*Gx1[72] + E1[26]*Gx1[86] + E1[30]*Gx1[100] + E1[34]*Gx1[114] + E1[38]*Gx1[128] + E1[42]*Gx1[142] + E1[46]*Gx1[156] + E1[50]*Gx1[170] + E1[54]*Gx1[184];
H101[31] += + E1[2]*Gx1[3] + E1[6]*Gx1[17] + E1[10]*Gx1[31] + E1[14]*Gx1[45] + E1[18]*Gx1[59] + E1[22]*Gx1[73] + E1[26]*Gx1[87] + E1[30]*Gx1[101] + E1[34]*Gx1[115] + E1[38]*Gx1[129] + E1[42]*Gx1[143] + E1[46]*Gx1[157] + E1[50]*Gx1[171] + E1[54]*Gx1[185];
H101[32] += + E1[2]*Gx1[4] + E1[6]*Gx1[18] + E1[10]*Gx1[32] + E1[14]*Gx1[46] + E1[18]*Gx1[60] + E1[22]*Gx1[74] + E1[26]*Gx1[88] + E1[30]*Gx1[102] + E1[34]*Gx1[116] + E1[38]*Gx1[130] + E1[42]*Gx1[144] + E1[46]*Gx1[158] + E1[50]*Gx1[172] + E1[54]*Gx1[186];
H101[33] += + E1[2]*Gx1[5] + E1[6]*Gx1[19] + E1[10]*Gx1[33] + E1[14]*Gx1[47] + E1[18]*Gx1[61] + E1[22]*Gx1[75] + E1[26]*Gx1[89] + E1[30]*Gx1[103] + E1[34]*Gx1[117] + E1[38]*Gx1[131] + E1[42]*Gx1[145] + E1[46]*Gx1[159] + E1[50]*Gx1[173] + E1[54]*Gx1[187];
H101[34] += + E1[2]*Gx1[6] + E1[6]*Gx1[20] + E1[10]*Gx1[34] + E1[14]*Gx1[48] + E1[18]*Gx1[62] + E1[22]*Gx1[76] + E1[26]*Gx1[90] + E1[30]*Gx1[104] + E1[34]*Gx1[118] + E1[38]*Gx1[132] + E1[42]*Gx1[146] + E1[46]*Gx1[160] + E1[50]*Gx1[174] + E1[54]*Gx1[188];
H101[35] += + E1[2]*Gx1[7] + E1[6]*Gx1[21] + E1[10]*Gx1[35] + E1[14]*Gx1[49] + E1[18]*Gx1[63] + E1[22]*Gx1[77] + E1[26]*Gx1[91] + E1[30]*Gx1[105] + E1[34]*Gx1[119] + E1[38]*Gx1[133] + E1[42]*Gx1[147] + E1[46]*Gx1[161] + E1[50]*Gx1[175] + E1[54]*Gx1[189];
H101[36] += + E1[2]*Gx1[8] + E1[6]*Gx1[22] + E1[10]*Gx1[36] + E1[14]*Gx1[50] + E1[18]*Gx1[64] + E1[22]*Gx1[78] + E1[26]*Gx1[92] + E1[30]*Gx1[106] + E1[34]*Gx1[120] + E1[38]*Gx1[134] + E1[42]*Gx1[148] + E1[46]*Gx1[162] + E1[50]*Gx1[176] + E1[54]*Gx1[190];
H101[37] += + E1[2]*Gx1[9] + E1[6]*Gx1[23] + E1[10]*Gx1[37] + E1[14]*Gx1[51] + E1[18]*Gx1[65] + E1[22]*Gx1[79] + E1[26]*Gx1[93] + E1[30]*Gx1[107] + E1[34]*Gx1[121] + E1[38]*Gx1[135] + E1[42]*Gx1[149] + E1[46]*Gx1[163] + E1[50]*Gx1[177] + E1[54]*Gx1[191];
H101[38] += + E1[2]*Gx1[10] + E1[6]*Gx1[24] + E1[10]*Gx1[38] + E1[14]*Gx1[52] + E1[18]*Gx1[66] + E1[22]*Gx1[80] + E1[26]*Gx1[94] + E1[30]*Gx1[108] + E1[34]*Gx1[122] + E1[38]*Gx1[136] + E1[42]*Gx1[150] + E1[46]*Gx1[164] + E1[50]*Gx1[178] + E1[54]*Gx1[192];
H101[39] += + E1[2]*Gx1[11] + E1[6]*Gx1[25] + E1[10]*Gx1[39] + E1[14]*Gx1[53] + E1[18]*Gx1[67] + E1[22]*Gx1[81] + E1[26]*Gx1[95] + E1[30]*Gx1[109] + E1[34]*Gx1[123] + E1[38]*Gx1[137] + E1[42]*Gx1[151] + E1[46]*Gx1[165] + E1[50]*Gx1[179] + E1[54]*Gx1[193];
H101[40] += + E1[2]*Gx1[12] + E1[6]*Gx1[26] + E1[10]*Gx1[40] + E1[14]*Gx1[54] + E1[18]*Gx1[68] + E1[22]*Gx1[82] + E1[26]*Gx1[96] + E1[30]*Gx1[110] + E1[34]*Gx1[124] + E1[38]*Gx1[138] + E1[42]*Gx1[152] + E1[46]*Gx1[166] + E1[50]*Gx1[180] + E1[54]*Gx1[194];
H101[41] += + E1[2]*Gx1[13] + E1[6]*Gx1[27] + E1[10]*Gx1[41] + E1[14]*Gx1[55] + E1[18]*Gx1[69] + E1[22]*Gx1[83] + E1[26]*Gx1[97] + E1[30]*Gx1[111] + E1[34]*Gx1[125] + E1[38]*Gx1[139] + E1[42]*Gx1[153] + E1[46]*Gx1[167] + E1[50]*Gx1[181] + E1[54]*Gx1[195];
H101[42] += + E1[3]*Gx1[0] + E1[7]*Gx1[14] + E1[11]*Gx1[28] + E1[15]*Gx1[42] + E1[19]*Gx1[56] + E1[23]*Gx1[70] + E1[27]*Gx1[84] + E1[31]*Gx1[98] + E1[35]*Gx1[112] + E1[39]*Gx1[126] + E1[43]*Gx1[140] + E1[47]*Gx1[154] + E1[51]*Gx1[168] + E1[55]*Gx1[182];
H101[43] += + E1[3]*Gx1[1] + E1[7]*Gx1[15] + E1[11]*Gx1[29] + E1[15]*Gx1[43] + E1[19]*Gx1[57] + E1[23]*Gx1[71] + E1[27]*Gx1[85] + E1[31]*Gx1[99] + E1[35]*Gx1[113] + E1[39]*Gx1[127] + E1[43]*Gx1[141] + E1[47]*Gx1[155] + E1[51]*Gx1[169] + E1[55]*Gx1[183];
H101[44] += + E1[3]*Gx1[2] + E1[7]*Gx1[16] + E1[11]*Gx1[30] + E1[15]*Gx1[44] + E1[19]*Gx1[58] + E1[23]*Gx1[72] + E1[27]*Gx1[86] + E1[31]*Gx1[100] + E1[35]*Gx1[114] + E1[39]*Gx1[128] + E1[43]*Gx1[142] + E1[47]*Gx1[156] + E1[51]*Gx1[170] + E1[55]*Gx1[184];
H101[45] += + E1[3]*Gx1[3] + E1[7]*Gx1[17] + E1[11]*Gx1[31] + E1[15]*Gx1[45] + E1[19]*Gx1[59] + E1[23]*Gx1[73] + E1[27]*Gx1[87] + E1[31]*Gx1[101] + E1[35]*Gx1[115] + E1[39]*Gx1[129] + E1[43]*Gx1[143] + E1[47]*Gx1[157] + E1[51]*Gx1[171] + E1[55]*Gx1[185];
H101[46] += + E1[3]*Gx1[4] + E1[7]*Gx1[18] + E1[11]*Gx1[32] + E1[15]*Gx1[46] + E1[19]*Gx1[60] + E1[23]*Gx1[74] + E1[27]*Gx1[88] + E1[31]*Gx1[102] + E1[35]*Gx1[116] + E1[39]*Gx1[130] + E1[43]*Gx1[144] + E1[47]*Gx1[158] + E1[51]*Gx1[172] + E1[55]*Gx1[186];
H101[47] += + E1[3]*Gx1[5] + E1[7]*Gx1[19] + E1[11]*Gx1[33] + E1[15]*Gx1[47] + E1[19]*Gx1[61] + E1[23]*Gx1[75] + E1[27]*Gx1[89] + E1[31]*Gx1[103] + E1[35]*Gx1[117] + E1[39]*Gx1[131] + E1[43]*Gx1[145] + E1[47]*Gx1[159] + E1[51]*Gx1[173] + E1[55]*Gx1[187];
H101[48] += + E1[3]*Gx1[6] + E1[7]*Gx1[20] + E1[11]*Gx1[34] + E1[15]*Gx1[48] + E1[19]*Gx1[62] + E1[23]*Gx1[76] + E1[27]*Gx1[90] + E1[31]*Gx1[104] + E1[35]*Gx1[118] + E1[39]*Gx1[132] + E1[43]*Gx1[146] + E1[47]*Gx1[160] + E1[51]*Gx1[174] + E1[55]*Gx1[188];
H101[49] += + E1[3]*Gx1[7] + E1[7]*Gx1[21] + E1[11]*Gx1[35] + E1[15]*Gx1[49] + E1[19]*Gx1[63] + E1[23]*Gx1[77] + E1[27]*Gx1[91] + E1[31]*Gx1[105] + E1[35]*Gx1[119] + E1[39]*Gx1[133] + E1[43]*Gx1[147] + E1[47]*Gx1[161] + E1[51]*Gx1[175] + E1[55]*Gx1[189];
H101[50] += + E1[3]*Gx1[8] + E1[7]*Gx1[22] + E1[11]*Gx1[36] + E1[15]*Gx1[50] + E1[19]*Gx1[64] + E1[23]*Gx1[78] + E1[27]*Gx1[92] + E1[31]*Gx1[106] + E1[35]*Gx1[120] + E1[39]*Gx1[134] + E1[43]*Gx1[148] + E1[47]*Gx1[162] + E1[51]*Gx1[176] + E1[55]*Gx1[190];
H101[51] += + E1[3]*Gx1[9] + E1[7]*Gx1[23] + E1[11]*Gx1[37] + E1[15]*Gx1[51] + E1[19]*Gx1[65] + E1[23]*Gx1[79] + E1[27]*Gx1[93] + E1[31]*Gx1[107] + E1[35]*Gx1[121] + E1[39]*Gx1[135] + E1[43]*Gx1[149] + E1[47]*Gx1[163] + E1[51]*Gx1[177] + E1[55]*Gx1[191];
H101[52] += + E1[3]*Gx1[10] + E1[7]*Gx1[24] + E1[11]*Gx1[38] + E1[15]*Gx1[52] + E1[19]*Gx1[66] + E1[23]*Gx1[80] + E1[27]*Gx1[94] + E1[31]*Gx1[108] + E1[35]*Gx1[122] + E1[39]*Gx1[136] + E1[43]*Gx1[150] + E1[47]*Gx1[164] + E1[51]*Gx1[178] + E1[55]*Gx1[192];
H101[53] += + E1[3]*Gx1[11] + E1[7]*Gx1[25] + E1[11]*Gx1[39] + E1[15]*Gx1[53] + E1[19]*Gx1[67] + E1[23]*Gx1[81] + E1[27]*Gx1[95] + E1[31]*Gx1[109] + E1[35]*Gx1[123] + E1[39]*Gx1[137] + E1[43]*Gx1[151] + E1[47]*Gx1[165] + E1[51]*Gx1[179] + E1[55]*Gx1[193];
H101[54] += + E1[3]*Gx1[12] + E1[7]*Gx1[26] + E1[11]*Gx1[40] + E1[15]*Gx1[54] + E1[19]*Gx1[68] + E1[23]*Gx1[82] + E1[27]*Gx1[96] + E1[31]*Gx1[110] + E1[35]*Gx1[124] + E1[39]*Gx1[138] + E1[43]*Gx1[152] + E1[47]*Gx1[166] + E1[51]*Gx1[180] + E1[55]*Gx1[194];
H101[55] += + E1[3]*Gx1[13] + E1[7]*Gx1[27] + E1[11]*Gx1[41] + E1[15]*Gx1[55] + E1[19]*Gx1[69] + E1[23]*Gx1[83] + E1[27]*Gx1[97] + E1[31]*Gx1[111] + E1[35]*Gx1[125] + E1[39]*Gx1[139] + E1[43]*Gx1[153] + E1[47]*Gx1[167] + E1[51]*Gx1[181] + E1[55]*Gx1[195];
}

void nmpc_zeroBlockH10( real_t* const H101 )
{
{ int lCopy; for (lCopy = 0; lCopy < 56; lCopy++) H101[ lCopy ] = 0; }
}

void nmpc_multEDu( real_t* const E1, real_t* const U1, real_t* const dNew )
{
dNew[0] += + E1[0]*U1[0] + E1[1]*U1[1] + E1[2]*U1[2] + E1[3]*U1[3];
dNew[1] += + E1[4]*U1[0] + E1[5]*U1[1] + E1[6]*U1[2] + E1[7]*U1[3];
dNew[2] += + E1[8]*U1[0] + E1[9]*U1[1] + E1[10]*U1[2] + E1[11]*U1[3];
dNew[3] += + E1[12]*U1[0] + E1[13]*U1[1] + E1[14]*U1[2] + E1[15]*U1[3];
dNew[4] += + E1[16]*U1[0] + E1[17]*U1[1] + E1[18]*U1[2] + E1[19]*U1[3];
dNew[5] += + E1[20]*U1[0] + E1[21]*U1[1] + E1[22]*U1[2] + E1[23]*U1[3];
dNew[6] += + E1[24]*U1[0] + E1[25]*U1[1] + E1[26]*U1[2] + E1[27]*U1[3];
dNew[7] += + E1[28]*U1[0] + E1[29]*U1[1] + E1[30]*U1[2] + E1[31]*U1[3];
dNew[8] += + E1[32]*U1[0] + E1[33]*U1[1] + E1[34]*U1[2] + E1[35]*U1[3];
dNew[9] += + E1[36]*U1[0] + E1[37]*U1[1] + E1[38]*U1[2] + E1[39]*U1[3];
dNew[10] += + E1[40]*U1[0] + E1[41]*U1[1] + E1[42]*U1[2] + E1[43]*U1[3];
dNew[11] += + E1[44]*U1[0] + E1[45]*U1[1] + E1[46]*U1[2] + E1[47]*U1[3];
dNew[12] += + E1[48]*U1[0] + E1[49]*U1[1] + E1[50]*U1[2] + E1[51]*U1[3];
dNew[13] += + E1[52]*U1[0] + E1[53]*U1[1] + E1[54]*U1[2] + E1[55]*U1[3];
}

void nmpc_zeroBlockH00(  )
{
nmpcWorkspace.H[0] = 0.0000000000000000e+00;
nmpcWorkspace.H[1] = 0.0000000000000000e+00;
nmpcWorkspace.H[2] = 0.0000000000000000e+00;
nmpcWorkspace.H[3] = 0.0000000000000000e+00;
nmpcWorkspace.H[4] = 0.0000000000000000e+00;
nmpcWorkspace.H[5] = 0.0000000000000000e+00;
nmpcWorkspace.H[6] = 0.0000000000000000e+00;
nmpcWorkspace.H[7] = 0.0000000000000000e+00;
nmpcWorkspace.H[8] = 0.0000000000000000e+00;
nmpcWorkspace.H[9] = 0.0000000000000000e+00;
nmpcWorkspace.H[10] = 0.0000000000000000e+00;
nmpcWorkspace.H[11] = 0.0000000000000000e+00;
nmpcWorkspace.H[12] = 0.0000000000000000e+00;
nmpcWorkspace.H[13] = 0.0000000000000000e+00;
nmpcWorkspace.H[214] = 0.0000000000000000e+00;
nmpcWorkspace.H[215] = 0.0000000000000000e+00;
nmpcWorkspace.H[216] = 0.0000000000000000e+00;
nmpcWorkspace.H[217] = 0.0000000000000000e+00;
nmpcWorkspace.H[218] = 0.0000000000000000e+00;
nmpcWorkspace.H[219] = 0.0000000000000000e+00;
nmpcWorkspace.H[220] = 0.0000000000000000e+00;
nmpcWorkspace.H[221] = 0.0000000000000000e+00;
nmpcWorkspace.H[222] = 0.0000000000000000e+00;
nmpcWorkspace.H[223] = 0.0000000000000000e+00;
nmpcWorkspace.H[224] = 0.0000000000000000e+00;
nmpcWorkspace.H[225] = 0.0000000000000000e+00;
nmpcWorkspace.H[226] = 0.0000000000000000e+00;
nmpcWorkspace.H[227] = 0.0000000000000000e+00;
nmpcWorkspace.H[428] = 0.0000000000000000e+00;
nmpcWorkspace.H[429] = 0.0000000000000000e+00;
nmpcWorkspace.H[430] = 0.0000000000000000e+00;
nmpcWorkspace.H[431] = 0.0000000000000000e+00;
nmpcWorkspace.H[432] = 0.0000000000000000e+00;
nmpcWorkspace.H[433] = 0.0000000000000000e+00;
nmpcWorkspace.H[434] = 0.0000000000000000e+00;
nmpcWorkspace.H[435] = 0.0000000000000000e+00;
nmpcWorkspace.H[436] = 0.0000000000000000e+00;
nmpcWorkspace.H[437] = 0.0000000000000000e+00;
nmpcWorkspace.H[438] = 0.0000000000000000e+00;
nmpcWorkspace.H[439] = 0.0000000000000000e+00;
nmpcWorkspace.H[440] = 0.0000000000000000e+00;
nmpcWorkspace.H[441] = 0.0000000000000000e+00;
nmpcWorkspace.H[642] = 0.0000000000000000e+00;
nmpcWorkspace.H[643] = 0.0000000000000000e+00;
nmpcWorkspace.H[644] = 0.0000000000000000e+00;
nmpcWorkspace.H[645] = 0.0000000000000000e+00;
nmpcWorkspace.H[646] = 0.0000000000000000e+00;
nmpcWorkspace.H[647] = 0.0000000000000000e+00;
nmpcWorkspace.H[648] = 0.0000000000000000e+00;
nmpcWorkspace.H[649] = 0.0000000000000000e+00;
nmpcWorkspace.H[650] = 0.0000000000000000e+00;
nmpcWorkspace.H[651] = 0.0000000000000000e+00;
nmpcWorkspace.H[652] = 0.0000000000000000e+00;
nmpcWorkspace.H[653] = 0.0000000000000000e+00;
nmpcWorkspace.H[654] = 0.0000000000000000e+00;
nmpcWorkspace.H[655] = 0.0000000000000000e+00;
nmpcWorkspace.H[856] = 0.0000000000000000e+00;
nmpcWorkspace.H[857] = 0.0000000000000000e+00;
nmpcWorkspace.H[858] = 0.0000000000000000e+00;
nmpcWorkspace.H[859] = 0.0000000000000000e+00;
nmpcWorkspace.H[860] = 0.0000000000000000e+00;
nmpcWorkspace.H[861] = 0.0000000000000000e+00;
nmpcWorkspace.H[862] = 0.0000000000000000e+00;
nmpcWorkspace.H[863] = 0.0000000000000000e+00;
nmpcWorkspace.H[864] = 0.0000000000000000e+00;
nmpcWorkspace.H[865] = 0.0000000000000000e+00;
nmpcWorkspace.H[866] = 0.0000000000000000e+00;
nmpcWorkspace.H[867] = 0.0000000000000000e+00;
nmpcWorkspace.H[868] = 0.0000000000000000e+00;
nmpcWorkspace.H[869] = 0.0000000000000000e+00;
nmpcWorkspace.H[1070] = 0.0000000000000000e+00;
nmpcWorkspace.H[1071] = 0.0000000000000000e+00;
nmpcWorkspace.H[1072] = 0.0000000000000000e+00;
nmpcWorkspace.H[1073] = 0.0000000000000000e+00;
nmpcWorkspace.H[1074] = 0.0000000000000000e+00;
nmpcWorkspace.H[1075] = 0.0000000000000000e+00;
nmpcWorkspace.H[1076] = 0.0000000000000000e+00;
nmpcWorkspace.H[1077] = 0.0000000000000000e+00;
nmpcWorkspace.H[1078] = 0.0000000000000000e+00;
nmpcWorkspace.H[1079] = 0.0000000000000000e+00;
nmpcWorkspace.H[1080] = 0.0000000000000000e+00;
nmpcWorkspace.H[1081] = 0.0000000000000000e+00;
nmpcWorkspace.H[1082] = 0.0000000000000000e+00;
nmpcWorkspace.H[1083] = 0.0000000000000000e+00;
nmpcWorkspace.H[1284] = 0.0000000000000000e+00;
nmpcWorkspace.H[1285] = 0.0000000000000000e+00;
nmpcWorkspace.H[1286] = 0.0000000000000000e+00;
nmpcWorkspace.H[1287] = 0.0000000000000000e+00;
nmpcWorkspace.H[1288] = 0.0000000000000000e+00;
nmpcWorkspace.H[1289] = 0.0000000000000000e+00;
nmpcWorkspace.H[1290] = 0.0000000000000000e+00;
nmpcWorkspace.H[1291] = 0.0000000000000000e+00;
nmpcWorkspace.H[1292] = 0.0000000000000000e+00;
nmpcWorkspace.H[1293] = 0.0000000000000000e+00;
nmpcWorkspace.H[1294] = 0.0000000000000000e+00;
nmpcWorkspace.H[1295] = 0.0000000000000000e+00;
nmpcWorkspace.H[1296] = 0.0000000000000000e+00;
nmpcWorkspace.H[1297] = 0.0000000000000000e+00;
nmpcWorkspace.H[1498] = 0.0000000000000000e+00;
nmpcWorkspace.H[1499] = 0.0000000000000000e+00;
nmpcWorkspace.H[1500] = 0.0000000000000000e+00;
nmpcWorkspace.H[1501] = 0.0000000000000000e+00;
nmpcWorkspace.H[1502] = 0.0000000000000000e+00;
nmpcWorkspace.H[1503] = 0.0000000000000000e+00;
nmpcWorkspace.H[1504] = 0.0000000000000000e+00;
nmpcWorkspace.H[1505] = 0.0000000000000000e+00;
nmpcWorkspace.H[1506] = 0.0000000000000000e+00;
nmpcWorkspace.H[1507] = 0.0000000000000000e+00;
nmpcWorkspace.H[1508] = 0.0000000000000000e+00;
nmpcWorkspace.H[1509] = 0.0000000000000000e+00;
nmpcWorkspace.H[1510] = 0.0000000000000000e+00;
nmpcWorkspace.H[1511] = 0.0000000000000000e+00;
nmpcWorkspace.H[1712] = 0.0000000000000000e+00;
nmpcWorkspace.H[1713] = 0.0000000000000000e+00;
nmpcWorkspace.H[1714] = 0.0000000000000000e+00;
nmpcWorkspace.H[1715] = 0.0000000000000000e+00;
nmpcWorkspace.H[1716] = 0.0000000000000000e+00;
nmpcWorkspace.H[1717] = 0.0000000000000000e+00;
nmpcWorkspace.H[1718] = 0.0000000000000000e+00;
nmpcWorkspace.H[1719] = 0.0000000000000000e+00;
nmpcWorkspace.H[1720] = 0.0000000000000000e+00;
nmpcWorkspace.H[1721] = 0.0000000000000000e+00;
nmpcWorkspace.H[1722] = 0.0000000000000000e+00;
nmpcWorkspace.H[1723] = 0.0000000000000000e+00;
nmpcWorkspace.H[1724] = 0.0000000000000000e+00;
nmpcWorkspace.H[1725] = 0.0000000000000000e+00;
nmpcWorkspace.H[1926] = 0.0000000000000000e+00;
nmpcWorkspace.H[1927] = 0.0000000000000000e+00;
nmpcWorkspace.H[1928] = 0.0000000000000000e+00;
nmpcWorkspace.H[1929] = 0.0000000000000000e+00;
nmpcWorkspace.H[1930] = 0.0000000000000000e+00;
nmpcWorkspace.H[1931] = 0.0000000000000000e+00;
nmpcWorkspace.H[1932] = 0.0000000000000000e+00;
nmpcWorkspace.H[1933] = 0.0000000000000000e+00;
nmpcWorkspace.H[1934] = 0.0000000000000000e+00;
nmpcWorkspace.H[1935] = 0.0000000000000000e+00;
nmpcWorkspace.H[1936] = 0.0000000000000000e+00;
nmpcWorkspace.H[1937] = 0.0000000000000000e+00;
nmpcWorkspace.H[1938] = 0.0000000000000000e+00;
nmpcWorkspace.H[1939] = 0.0000000000000000e+00;
nmpcWorkspace.H[2140] = 0.0000000000000000e+00;
nmpcWorkspace.H[2141] = 0.0000000000000000e+00;
nmpcWorkspace.H[2142] = 0.0000000000000000e+00;
nmpcWorkspace.H[2143] = 0.0000000000000000e+00;
nmpcWorkspace.H[2144] = 0.0000000000000000e+00;
nmpcWorkspace.H[2145] = 0.0000000000000000e+00;
nmpcWorkspace.H[2146] = 0.0000000000000000e+00;
nmpcWorkspace.H[2147] = 0.0000000000000000e+00;
nmpcWorkspace.H[2148] = 0.0000000000000000e+00;
nmpcWorkspace.H[2149] = 0.0000000000000000e+00;
nmpcWorkspace.H[2150] = 0.0000000000000000e+00;
nmpcWorkspace.H[2151] = 0.0000000000000000e+00;
nmpcWorkspace.H[2152] = 0.0000000000000000e+00;
nmpcWorkspace.H[2153] = 0.0000000000000000e+00;
nmpcWorkspace.H[2354] = 0.0000000000000000e+00;
nmpcWorkspace.H[2355] = 0.0000000000000000e+00;
nmpcWorkspace.H[2356] = 0.0000000000000000e+00;
nmpcWorkspace.H[2357] = 0.0000000000000000e+00;
nmpcWorkspace.H[2358] = 0.0000000000000000e+00;
nmpcWorkspace.H[2359] = 0.0000000000000000e+00;
nmpcWorkspace.H[2360] = 0.0000000000000000e+00;
nmpcWorkspace.H[2361] = 0.0000000000000000e+00;
nmpcWorkspace.H[2362] = 0.0000000000000000e+00;
nmpcWorkspace.H[2363] = 0.0000000000000000e+00;
nmpcWorkspace.H[2364] = 0.0000000000000000e+00;
nmpcWorkspace.H[2365] = 0.0000000000000000e+00;
nmpcWorkspace.H[2366] = 0.0000000000000000e+00;
nmpcWorkspace.H[2367] = 0.0000000000000000e+00;
nmpcWorkspace.H[2568] = 0.0000000000000000e+00;
nmpcWorkspace.H[2569] = 0.0000000000000000e+00;
nmpcWorkspace.H[2570] = 0.0000000000000000e+00;
nmpcWorkspace.H[2571] = 0.0000000000000000e+00;
nmpcWorkspace.H[2572] = 0.0000000000000000e+00;
nmpcWorkspace.H[2573] = 0.0000000000000000e+00;
nmpcWorkspace.H[2574] = 0.0000000000000000e+00;
nmpcWorkspace.H[2575] = 0.0000000000000000e+00;
nmpcWorkspace.H[2576] = 0.0000000000000000e+00;
nmpcWorkspace.H[2577] = 0.0000000000000000e+00;
nmpcWorkspace.H[2578] = 0.0000000000000000e+00;
nmpcWorkspace.H[2579] = 0.0000000000000000e+00;
nmpcWorkspace.H[2580] = 0.0000000000000000e+00;
nmpcWorkspace.H[2581] = 0.0000000000000000e+00;
nmpcWorkspace.H[2782] = 0.0000000000000000e+00;
nmpcWorkspace.H[2783] = 0.0000000000000000e+00;
nmpcWorkspace.H[2784] = 0.0000000000000000e+00;
nmpcWorkspace.H[2785] = 0.0000000000000000e+00;
nmpcWorkspace.H[2786] = 0.0000000000000000e+00;
nmpcWorkspace.H[2787] = 0.0000000000000000e+00;
nmpcWorkspace.H[2788] = 0.0000000000000000e+00;
nmpcWorkspace.H[2789] = 0.0000000000000000e+00;
nmpcWorkspace.H[2790] = 0.0000000000000000e+00;
nmpcWorkspace.H[2791] = 0.0000000000000000e+00;
nmpcWorkspace.H[2792] = 0.0000000000000000e+00;
nmpcWorkspace.H[2793] = 0.0000000000000000e+00;
nmpcWorkspace.H[2794] = 0.0000000000000000e+00;
nmpcWorkspace.H[2795] = 0.0000000000000000e+00;
}

void nmpc_multCTQC( real_t* const Gx1, real_t* const Gx2 )
{
nmpcWorkspace.H[0] += + Gx1[0]*Gx2[0] + Gx1[14]*Gx2[14] + Gx1[28]*Gx2[28] + Gx1[42]*Gx2[42] + Gx1[56]*Gx2[56] + Gx1[70]*Gx2[70] + Gx1[84]*Gx2[84] + Gx1[98]*Gx2[98] + Gx1[112]*Gx2[112] + Gx1[126]*Gx2[126] + Gx1[140]*Gx2[140] + Gx1[154]*Gx2[154] + Gx1[168]*Gx2[168] + Gx1[182]*Gx2[182];
nmpcWorkspace.H[1] += + Gx1[0]*Gx2[1] + Gx1[14]*Gx2[15] + Gx1[28]*Gx2[29] + Gx1[42]*Gx2[43] + Gx1[56]*Gx2[57] + Gx1[70]*Gx2[71] + Gx1[84]*Gx2[85] + Gx1[98]*Gx2[99] + Gx1[112]*Gx2[113] + Gx1[126]*Gx2[127] + Gx1[140]*Gx2[141] + Gx1[154]*Gx2[155] + Gx1[168]*Gx2[169] + Gx1[182]*Gx2[183];
nmpcWorkspace.H[2] += + Gx1[0]*Gx2[2] + Gx1[14]*Gx2[16] + Gx1[28]*Gx2[30] + Gx1[42]*Gx2[44] + Gx1[56]*Gx2[58] + Gx1[70]*Gx2[72] + Gx1[84]*Gx2[86] + Gx1[98]*Gx2[100] + Gx1[112]*Gx2[114] + Gx1[126]*Gx2[128] + Gx1[140]*Gx2[142] + Gx1[154]*Gx2[156] + Gx1[168]*Gx2[170] + Gx1[182]*Gx2[184];
nmpcWorkspace.H[3] += + Gx1[0]*Gx2[3] + Gx1[14]*Gx2[17] + Gx1[28]*Gx2[31] + Gx1[42]*Gx2[45] + Gx1[56]*Gx2[59] + Gx1[70]*Gx2[73] + Gx1[84]*Gx2[87] + Gx1[98]*Gx2[101] + Gx1[112]*Gx2[115] + Gx1[126]*Gx2[129] + Gx1[140]*Gx2[143] + Gx1[154]*Gx2[157] + Gx1[168]*Gx2[171] + Gx1[182]*Gx2[185];
nmpcWorkspace.H[4] += + Gx1[0]*Gx2[4] + Gx1[14]*Gx2[18] + Gx1[28]*Gx2[32] + Gx1[42]*Gx2[46] + Gx1[56]*Gx2[60] + Gx1[70]*Gx2[74] + Gx1[84]*Gx2[88] + Gx1[98]*Gx2[102] + Gx1[112]*Gx2[116] + Gx1[126]*Gx2[130] + Gx1[140]*Gx2[144] + Gx1[154]*Gx2[158] + Gx1[168]*Gx2[172] + Gx1[182]*Gx2[186];
nmpcWorkspace.H[5] += + Gx1[0]*Gx2[5] + Gx1[14]*Gx2[19] + Gx1[28]*Gx2[33] + Gx1[42]*Gx2[47] + Gx1[56]*Gx2[61] + Gx1[70]*Gx2[75] + Gx1[84]*Gx2[89] + Gx1[98]*Gx2[103] + Gx1[112]*Gx2[117] + Gx1[126]*Gx2[131] + Gx1[140]*Gx2[145] + Gx1[154]*Gx2[159] + Gx1[168]*Gx2[173] + Gx1[182]*Gx2[187];
nmpcWorkspace.H[6] += + Gx1[0]*Gx2[6] + Gx1[14]*Gx2[20] + Gx1[28]*Gx2[34] + Gx1[42]*Gx2[48] + Gx1[56]*Gx2[62] + Gx1[70]*Gx2[76] + Gx1[84]*Gx2[90] + Gx1[98]*Gx2[104] + Gx1[112]*Gx2[118] + Gx1[126]*Gx2[132] + Gx1[140]*Gx2[146] + Gx1[154]*Gx2[160] + Gx1[168]*Gx2[174] + Gx1[182]*Gx2[188];
nmpcWorkspace.H[7] += + Gx1[0]*Gx2[7] + Gx1[14]*Gx2[21] + Gx1[28]*Gx2[35] + Gx1[42]*Gx2[49] + Gx1[56]*Gx2[63] + Gx1[70]*Gx2[77] + Gx1[84]*Gx2[91] + Gx1[98]*Gx2[105] + Gx1[112]*Gx2[119] + Gx1[126]*Gx2[133] + Gx1[140]*Gx2[147] + Gx1[154]*Gx2[161] + Gx1[168]*Gx2[175] + Gx1[182]*Gx2[189];
nmpcWorkspace.H[8] += + Gx1[0]*Gx2[8] + Gx1[14]*Gx2[22] + Gx1[28]*Gx2[36] + Gx1[42]*Gx2[50] + Gx1[56]*Gx2[64] + Gx1[70]*Gx2[78] + Gx1[84]*Gx2[92] + Gx1[98]*Gx2[106] + Gx1[112]*Gx2[120] + Gx1[126]*Gx2[134] + Gx1[140]*Gx2[148] + Gx1[154]*Gx2[162] + Gx1[168]*Gx2[176] + Gx1[182]*Gx2[190];
nmpcWorkspace.H[9] += + Gx1[0]*Gx2[9] + Gx1[14]*Gx2[23] + Gx1[28]*Gx2[37] + Gx1[42]*Gx2[51] + Gx1[56]*Gx2[65] + Gx1[70]*Gx2[79] + Gx1[84]*Gx2[93] + Gx1[98]*Gx2[107] + Gx1[112]*Gx2[121] + Gx1[126]*Gx2[135] + Gx1[140]*Gx2[149] + Gx1[154]*Gx2[163] + Gx1[168]*Gx2[177] + Gx1[182]*Gx2[191];
nmpcWorkspace.H[10] += + Gx1[0]*Gx2[10] + Gx1[14]*Gx2[24] + Gx1[28]*Gx2[38] + Gx1[42]*Gx2[52] + Gx1[56]*Gx2[66] + Gx1[70]*Gx2[80] + Gx1[84]*Gx2[94] + Gx1[98]*Gx2[108] + Gx1[112]*Gx2[122] + Gx1[126]*Gx2[136] + Gx1[140]*Gx2[150] + Gx1[154]*Gx2[164] + Gx1[168]*Gx2[178] + Gx1[182]*Gx2[192];
nmpcWorkspace.H[11] += + Gx1[0]*Gx2[11] + Gx1[14]*Gx2[25] + Gx1[28]*Gx2[39] + Gx1[42]*Gx2[53] + Gx1[56]*Gx2[67] + Gx1[70]*Gx2[81] + Gx1[84]*Gx2[95] + Gx1[98]*Gx2[109] + Gx1[112]*Gx2[123] + Gx1[126]*Gx2[137] + Gx1[140]*Gx2[151] + Gx1[154]*Gx2[165] + Gx1[168]*Gx2[179] + Gx1[182]*Gx2[193];
nmpcWorkspace.H[12] += + Gx1[0]*Gx2[12] + Gx1[14]*Gx2[26] + Gx1[28]*Gx2[40] + Gx1[42]*Gx2[54] + Gx1[56]*Gx2[68] + Gx1[70]*Gx2[82] + Gx1[84]*Gx2[96] + Gx1[98]*Gx2[110] + Gx1[112]*Gx2[124] + Gx1[126]*Gx2[138] + Gx1[140]*Gx2[152] + Gx1[154]*Gx2[166] + Gx1[168]*Gx2[180] + Gx1[182]*Gx2[194];
nmpcWorkspace.H[13] += + Gx1[0]*Gx2[13] + Gx1[14]*Gx2[27] + Gx1[28]*Gx2[41] + Gx1[42]*Gx2[55] + Gx1[56]*Gx2[69] + Gx1[70]*Gx2[83] + Gx1[84]*Gx2[97] + Gx1[98]*Gx2[111] + Gx1[112]*Gx2[125] + Gx1[126]*Gx2[139] + Gx1[140]*Gx2[153] + Gx1[154]*Gx2[167] + Gx1[168]*Gx2[181] + Gx1[182]*Gx2[195];
nmpcWorkspace.H[214] += + Gx1[1]*Gx2[0] + Gx1[15]*Gx2[14] + Gx1[29]*Gx2[28] + Gx1[43]*Gx2[42] + Gx1[57]*Gx2[56] + Gx1[71]*Gx2[70] + Gx1[85]*Gx2[84] + Gx1[99]*Gx2[98] + Gx1[113]*Gx2[112] + Gx1[127]*Gx2[126] + Gx1[141]*Gx2[140] + Gx1[155]*Gx2[154] + Gx1[169]*Gx2[168] + Gx1[183]*Gx2[182];
nmpcWorkspace.H[215] += + Gx1[1]*Gx2[1] + Gx1[15]*Gx2[15] + Gx1[29]*Gx2[29] + Gx1[43]*Gx2[43] + Gx1[57]*Gx2[57] + Gx1[71]*Gx2[71] + Gx1[85]*Gx2[85] + Gx1[99]*Gx2[99] + Gx1[113]*Gx2[113] + Gx1[127]*Gx2[127] + Gx1[141]*Gx2[141] + Gx1[155]*Gx2[155] + Gx1[169]*Gx2[169] + Gx1[183]*Gx2[183];
nmpcWorkspace.H[216] += + Gx1[1]*Gx2[2] + Gx1[15]*Gx2[16] + Gx1[29]*Gx2[30] + Gx1[43]*Gx2[44] + Gx1[57]*Gx2[58] + Gx1[71]*Gx2[72] + Gx1[85]*Gx2[86] + Gx1[99]*Gx2[100] + Gx1[113]*Gx2[114] + Gx1[127]*Gx2[128] + Gx1[141]*Gx2[142] + Gx1[155]*Gx2[156] + Gx1[169]*Gx2[170] + Gx1[183]*Gx2[184];
nmpcWorkspace.H[217] += + Gx1[1]*Gx2[3] + Gx1[15]*Gx2[17] + Gx1[29]*Gx2[31] + Gx1[43]*Gx2[45] + Gx1[57]*Gx2[59] + Gx1[71]*Gx2[73] + Gx1[85]*Gx2[87] + Gx1[99]*Gx2[101] + Gx1[113]*Gx2[115] + Gx1[127]*Gx2[129] + Gx1[141]*Gx2[143] + Gx1[155]*Gx2[157] + Gx1[169]*Gx2[171] + Gx1[183]*Gx2[185];
nmpcWorkspace.H[218] += + Gx1[1]*Gx2[4] + Gx1[15]*Gx2[18] + Gx1[29]*Gx2[32] + Gx1[43]*Gx2[46] + Gx1[57]*Gx2[60] + Gx1[71]*Gx2[74] + Gx1[85]*Gx2[88] + Gx1[99]*Gx2[102] + Gx1[113]*Gx2[116] + Gx1[127]*Gx2[130] + Gx1[141]*Gx2[144] + Gx1[155]*Gx2[158] + Gx1[169]*Gx2[172] + Gx1[183]*Gx2[186];
nmpcWorkspace.H[219] += + Gx1[1]*Gx2[5] + Gx1[15]*Gx2[19] + Gx1[29]*Gx2[33] + Gx1[43]*Gx2[47] + Gx1[57]*Gx2[61] + Gx1[71]*Gx2[75] + Gx1[85]*Gx2[89] + Gx1[99]*Gx2[103] + Gx1[113]*Gx2[117] + Gx1[127]*Gx2[131] + Gx1[141]*Gx2[145] + Gx1[155]*Gx2[159] + Gx1[169]*Gx2[173] + Gx1[183]*Gx2[187];
nmpcWorkspace.H[220] += + Gx1[1]*Gx2[6] + Gx1[15]*Gx2[20] + Gx1[29]*Gx2[34] + Gx1[43]*Gx2[48] + Gx1[57]*Gx2[62] + Gx1[71]*Gx2[76] + Gx1[85]*Gx2[90] + Gx1[99]*Gx2[104] + Gx1[113]*Gx2[118] + Gx1[127]*Gx2[132] + Gx1[141]*Gx2[146] + Gx1[155]*Gx2[160] + Gx1[169]*Gx2[174] + Gx1[183]*Gx2[188];
nmpcWorkspace.H[221] += + Gx1[1]*Gx2[7] + Gx1[15]*Gx2[21] + Gx1[29]*Gx2[35] + Gx1[43]*Gx2[49] + Gx1[57]*Gx2[63] + Gx1[71]*Gx2[77] + Gx1[85]*Gx2[91] + Gx1[99]*Gx2[105] + Gx1[113]*Gx2[119] + Gx1[127]*Gx2[133] + Gx1[141]*Gx2[147] + Gx1[155]*Gx2[161] + Gx1[169]*Gx2[175] + Gx1[183]*Gx2[189];
nmpcWorkspace.H[222] += + Gx1[1]*Gx2[8] + Gx1[15]*Gx2[22] + Gx1[29]*Gx2[36] + Gx1[43]*Gx2[50] + Gx1[57]*Gx2[64] + Gx1[71]*Gx2[78] + Gx1[85]*Gx2[92] + Gx1[99]*Gx2[106] + Gx1[113]*Gx2[120] + Gx1[127]*Gx2[134] + Gx1[141]*Gx2[148] + Gx1[155]*Gx2[162] + Gx1[169]*Gx2[176] + Gx1[183]*Gx2[190];
nmpcWorkspace.H[223] += + Gx1[1]*Gx2[9] + Gx1[15]*Gx2[23] + Gx1[29]*Gx2[37] + Gx1[43]*Gx2[51] + Gx1[57]*Gx2[65] + Gx1[71]*Gx2[79] + Gx1[85]*Gx2[93] + Gx1[99]*Gx2[107] + Gx1[113]*Gx2[121] + Gx1[127]*Gx2[135] + Gx1[141]*Gx2[149] + Gx1[155]*Gx2[163] + Gx1[169]*Gx2[177] + Gx1[183]*Gx2[191];
nmpcWorkspace.H[224] += + Gx1[1]*Gx2[10] + Gx1[15]*Gx2[24] + Gx1[29]*Gx2[38] + Gx1[43]*Gx2[52] + Gx1[57]*Gx2[66] + Gx1[71]*Gx2[80] + Gx1[85]*Gx2[94] + Gx1[99]*Gx2[108] + Gx1[113]*Gx2[122] + Gx1[127]*Gx2[136] + Gx1[141]*Gx2[150] + Gx1[155]*Gx2[164] + Gx1[169]*Gx2[178] + Gx1[183]*Gx2[192];
nmpcWorkspace.H[225] += + Gx1[1]*Gx2[11] + Gx1[15]*Gx2[25] + Gx1[29]*Gx2[39] + Gx1[43]*Gx2[53] + Gx1[57]*Gx2[67] + Gx1[71]*Gx2[81] + Gx1[85]*Gx2[95] + Gx1[99]*Gx2[109] + Gx1[113]*Gx2[123] + Gx1[127]*Gx2[137] + Gx1[141]*Gx2[151] + Gx1[155]*Gx2[165] + Gx1[169]*Gx2[179] + Gx1[183]*Gx2[193];
nmpcWorkspace.H[226] += + Gx1[1]*Gx2[12] + Gx1[15]*Gx2[26] + Gx1[29]*Gx2[40] + Gx1[43]*Gx2[54] + Gx1[57]*Gx2[68] + Gx1[71]*Gx2[82] + Gx1[85]*Gx2[96] + Gx1[99]*Gx2[110] + Gx1[113]*Gx2[124] + Gx1[127]*Gx2[138] + Gx1[141]*Gx2[152] + Gx1[155]*Gx2[166] + Gx1[169]*Gx2[180] + Gx1[183]*Gx2[194];
nmpcWorkspace.H[227] += + Gx1[1]*Gx2[13] + Gx1[15]*Gx2[27] + Gx1[29]*Gx2[41] + Gx1[43]*Gx2[55] + Gx1[57]*Gx2[69] + Gx1[71]*Gx2[83] + Gx1[85]*Gx2[97] + Gx1[99]*Gx2[111] + Gx1[113]*Gx2[125] + Gx1[127]*Gx2[139] + Gx1[141]*Gx2[153] + Gx1[155]*Gx2[167] + Gx1[169]*Gx2[181] + Gx1[183]*Gx2[195];
nmpcWorkspace.H[428] += + Gx1[2]*Gx2[0] + Gx1[16]*Gx2[14] + Gx1[30]*Gx2[28] + Gx1[44]*Gx2[42] + Gx1[58]*Gx2[56] + Gx1[72]*Gx2[70] + Gx1[86]*Gx2[84] + Gx1[100]*Gx2[98] + Gx1[114]*Gx2[112] + Gx1[128]*Gx2[126] + Gx1[142]*Gx2[140] + Gx1[156]*Gx2[154] + Gx1[170]*Gx2[168] + Gx1[184]*Gx2[182];
nmpcWorkspace.H[429] += + Gx1[2]*Gx2[1] + Gx1[16]*Gx2[15] + Gx1[30]*Gx2[29] + Gx1[44]*Gx2[43] + Gx1[58]*Gx2[57] + Gx1[72]*Gx2[71] + Gx1[86]*Gx2[85] + Gx1[100]*Gx2[99] + Gx1[114]*Gx2[113] + Gx1[128]*Gx2[127] + Gx1[142]*Gx2[141] + Gx1[156]*Gx2[155] + Gx1[170]*Gx2[169] + Gx1[184]*Gx2[183];
nmpcWorkspace.H[430] += + Gx1[2]*Gx2[2] + Gx1[16]*Gx2[16] + Gx1[30]*Gx2[30] + Gx1[44]*Gx2[44] + Gx1[58]*Gx2[58] + Gx1[72]*Gx2[72] + Gx1[86]*Gx2[86] + Gx1[100]*Gx2[100] + Gx1[114]*Gx2[114] + Gx1[128]*Gx2[128] + Gx1[142]*Gx2[142] + Gx1[156]*Gx2[156] + Gx1[170]*Gx2[170] + Gx1[184]*Gx2[184];
nmpcWorkspace.H[431] += + Gx1[2]*Gx2[3] + Gx1[16]*Gx2[17] + Gx1[30]*Gx2[31] + Gx1[44]*Gx2[45] + Gx1[58]*Gx2[59] + Gx1[72]*Gx2[73] + Gx1[86]*Gx2[87] + Gx1[100]*Gx2[101] + Gx1[114]*Gx2[115] + Gx1[128]*Gx2[129] + Gx1[142]*Gx2[143] + Gx1[156]*Gx2[157] + Gx1[170]*Gx2[171] + Gx1[184]*Gx2[185];
nmpcWorkspace.H[432] += + Gx1[2]*Gx2[4] + Gx1[16]*Gx2[18] + Gx1[30]*Gx2[32] + Gx1[44]*Gx2[46] + Gx1[58]*Gx2[60] + Gx1[72]*Gx2[74] + Gx1[86]*Gx2[88] + Gx1[100]*Gx2[102] + Gx1[114]*Gx2[116] + Gx1[128]*Gx2[130] + Gx1[142]*Gx2[144] + Gx1[156]*Gx2[158] + Gx1[170]*Gx2[172] + Gx1[184]*Gx2[186];
nmpcWorkspace.H[433] += + Gx1[2]*Gx2[5] + Gx1[16]*Gx2[19] + Gx1[30]*Gx2[33] + Gx1[44]*Gx2[47] + Gx1[58]*Gx2[61] + Gx1[72]*Gx2[75] + Gx1[86]*Gx2[89] + Gx1[100]*Gx2[103] + Gx1[114]*Gx2[117] + Gx1[128]*Gx2[131] + Gx1[142]*Gx2[145] + Gx1[156]*Gx2[159] + Gx1[170]*Gx2[173] + Gx1[184]*Gx2[187];
nmpcWorkspace.H[434] += + Gx1[2]*Gx2[6] + Gx1[16]*Gx2[20] + Gx1[30]*Gx2[34] + Gx1[44]*Gx2[48] + Gx1[58]*Gx2[62] + Gx1[72]*Gx2[76] + Gx1[86]*Gx2[90] + Gx1[100]*Gx2[104] + Gx1[114]*Gx2[118] + Gx1[128]*Gx2[132] + Gx1[142]*Gx2[146] + Gx1[156]*Gx2[160] + Gx1[170]*Gx2[174] + Gx1[184]*Gx2[188];
nmpcWorkspace.H[435] += + Gx1[2]*Gx2[7] + Gx1[16]*Gx2[21] + Gx1[30]*Gx2[35] + Gx1[44]*Gx2[49] + Gx1[58]*Gx2[63] + Gx1[72]*Gx2[77] + Gx1[86]*Gx2[91] + Gx1[100]*Gx2[105] + Gx1[114]*Gx2[119] + Gx1[128]*Gx2[133] + Gx1[142]*Gx2[147] + Gx1[156]*Gx2[161] + Gx1[170]*Gx2[175] + Gx1[184]*Gx2[189];
nmpcWorkspace.H[436] += + Gx1[2]*Gx2[8] + Gx1[16]*Gx2[22] + Gx1[30]*Gx2[36] + Gx1[44]*Gx2[50] + Gx1[58]*Gx2[64] + Gx1[72]*Gx2[78] + Gx1[86]*Gx2[92] + Gx1[100]*Gx2[106] + Gx1[114]*Gx2[120] + Gx1[128]*Gx2[134] + Gx1[142]*Gx2[148] + Gx1[156]*Gx2[162] + Gx1[170]*Gx2[176] + Gx1[184]*Gx2[190];
nmpcWorkspace.H[437] += + Gx1[2]*Gx2[9] + Gx1[16]*Gx2[23] + Gx1[30]*Gx2[37] + Gx1[44]*Gx2[51] + Gx1[58]*Gx2[65] + Gx1[72]*Gx2[79] + Gx1[86]*Gx2[93] + Gx1[100]*Gx2[107] + Gx1[114]*Gx2[121] + Gx1[128]*Gx2[135] + Gx1[142]*Gx2[149] + Gx1[156]*Gx2[163] + Gx1[170]*Gx2[177] + Gx1[184]*Gx2[191];
nmpcWorkspace.H[438] += + Gx1[2]*Gx2[10] + Gx1[16]*Gx2[24] + Gx1[30]*Gx2[38] + Gx1[44]*Gx2[52] + Gx1[58]*Gx2[66] + Gx1[72]*Gx2[80] + Gx1[86]*Gx2[94] + Gx1[100]*Gx2[108] + Gx1[114]*Gx2[122] + Gx1[128]*Gx2[136] + Gx1[142]*Gx2[150] + Gx1[156]*Gx2[164] + Gx1[170]*Gx2[178] + Gx1[184]*Gx2[192];
nmpcWorkspace.H[439] += + Gx1[2]*Gx2[11] + Gx1[16]*Gx2[25] + Gx1[30]*Gx2[39] + Gx1[44]*Gx2[53] + Gx1[58]*Gx2[67] + Gx1[72]*Gx2[81] + Gx1[86]*Gx2[95] + Gx1[100]*Gx2[109] + Gx1[114]*Gx2[123] + Gx1[128]*Gx2[137] + Gx1[142]*Gx2[151] + Gx1[156]*Gx2[165] + Gx1[170]*Gx2[179] + Gx1[184]*Gx2[193];
nmpcWorkspace.H[440] += + Gx1[2]*Gx2[12] + Gx1[16]*Gx2[26] + Gx1[30]*Gx2[40] + Gx1[44]*Gx2[54] + Gx1[58]*Gx2[68] + Gx1[72]*Gx2[82] + Gx1[86]*Gx2[96] + Gx1[100]*Gx2[110] + Gx1[114]*Gx2[124] + Gx1[128]*Gx2[138] + Gx1[142]*Gx2[152] + Gx1[156]*Gx2[166] + Gx1[170]*Gx2[180] + Gx1[184]*Gx2[194];
nmpcWorkspace.H[441] += + Gx1[2]*Gx2[13] + Gx1[16]*Gx2[27] + Gx1[30]*Gx2[41] + Gx1[44]*Gx2[55] + Gx1[58]*Gx2[69] + Gx1[72]*Gx2[83] + Gx1[86]*Gx2[97] + Gx1[100]*Gx2[111] + Gx1[114]*Gx2[125] + Gx1[128]*Gx2[139] + Gx1[142]*Gx2[153] + Gx1[156]*Gx2[167] + Gx1[170]*Gx2[181] + Gx1[184]*Gx2[195];
nmpcWorkspace.H[642] += + Gx1[3]*Gx2[0] + Gx1[17]*Gx2[14] + Gx1[31]*Gx2[28] + Gx1[45]*Gx2[42] + Gx1[59]*Gx2[56] + Gx1[73]*Gx2[70] + Gx1[87]*Gx2[84] + Gx1[101]*Gx2[98] + Gx1[115]*Gx2[112] + Gx1[129]*Gx2[126] + Gx1[143]*Gx2[140] + Gx1[157]*Gx2[154] + Gx1[171]*Gx2[168] + Gx1[185]*Gx2[182];
nmpcWorkspace.H[643] += + Gx1[3]*Gx2[1] + Gx1[17]*Gx2[15] + Gx1[31]*Gx2[29] + Gx1[45]*Gx2[43] + Gx1[59]*Gx2[57] + Gx1[73]*Gx2[71] + Gx1[87]*Gx2[85] + Gx1[101]*Gx2[99] + Gx1[115]*Gx2[113] + Gx1[129]*Gx2[127] + Gx1[143]*Gx2[141] + Gx1[157]*Gx2[155] + Gx1[171]*Gx2[169] + Gx1[185]*Gx2[183];
nmpcWorkspace.H[644] += + Gx1[3]*Gx2[2] + Gx1[17]*Gx2[16] + Gx1[31]*Gx2[30] + Gx1[45]*Gx2[44] + Gx1[59]*Gx2[58] + Gx1[73]*Gx2[72] + Gx1[87]*Gx2[86] + Gx1[101]*Gx2[100] + Gx1[115]*Gx2[114] + Gx1[129]*Gx2[128] + Gx1[143]*Gx2[142] + Gx1[157]*Gx2[156] + Gx1[171]*Gx2[170] + Gx1[185]*Gx2[184];
nmpcWorkspace.H[645] += + Gx1[3]*Gx2[3] + Gx1[17]*Gx2[17] + Gx1[31]*Gx2[31] + Gx1[45]*Gx2[45] + Gx1[59]*Gx2[59] + Gx1[73]*Gx2[73] + Gx1[87]*Gx2[87] + Gx1[101]*Gx2[101] + Gx1[115]*Gx2[115] + Gx1[129]*Gx2[129] + Gx1[143]*Gx2[143] + Gx1[157]*Gx2[157] + Gx1[171]*Gx2[171] + Gx1[185]*Gx2[185];
nmpcWorkspace.H[646] += + Gx1[3]*Gx2[4] + Gx1[17]*Gx2[18] + Gx1[31]*Gx2[32] + Gx1[45]*Gx2[46] + Gx1[59]*Gx2[60] + Gx1[73]*Gx2[74] + Gx1[87]*Gx2[88] + Gx1[101]*Gx2[102] + Gx1[115]*Gx2[116] + Gx1[129]*Gx2[130] + Gx1[143]*Gx2[144] + Gx1[157]*Gx2[158] + Gx1[171]*Gx2[172] + Gx1[185]*Gx2[186];
nmpcWorkspace.H[647] += + Gx1[3]*Gx2[5] + Gx1[17]*Gx2[19] + Gx1[31]*Gx2[33] + Gx1[45]*Gx2[47] + Gx1[59]*Gx2[61] + Gx1[73]*Gx2[75] + Gx1[87]*Gx2[89] + Gx1[101]*Gx2[103] + Gx1[115]*Gx2[117] + Gx1[129]*Gx2[131] + Gx1[143]*Gx2[145] + Gx1[157]*Gx2[159] + Gx1[171]*Gx2[173] + Gx1[185]*Gx2[187];
nmpcWorkspace.H[648] += + Gx1[3]*Gx2[6] + Gx1[17]*Gx2[20] + Gx1[31]*Gx2[34] + Gx1[45]*Gx2[48] + Gx1[59]*Gx2[62] + Gx1[73]*Gx2[76] + Gx1[87]*Gx2[90] + Gx1[101]*Gx2[104] + Gx1[115]*Gx2[118] + Gx1[129]*Gx2[132] + Gx1[143]*Gx2[146] + Gx1[157]*Gx2[160] + Gx1[171]*Gx2[174] + Gx1[185]*Gx2[188];
nmpcWorkspace.H[649] += + Gx1[3]*Gx2[7] + Gx1[17]*Gx2[21] + Gx1[31]*Gx2[35] + Gx1[45]*Gx2[49] + Gx1[59]*Gx2[63] + Gx1[73]*Gx2[77] + Gx1[87]*Gx2[91] + Gx1[101]*Gx2[105] + Gx1[115]*Gx2[119] + Gx1[129]*Gx2[133] + Gx1[143]*Gx2[147] + Gx1[157]*Gx2[161] + Gx1[171]*Gx2[175] + Gx1[185]*Gx2[189];
nmpcWorkspace.H[650] += + Gx1[3]*Gx2[8] + Gx1[17]*Gx2[22] + Gx1[31]*Gx2[36] + Gx1[45]*Gx2[50] + Gx1[59]*Gx2[64] + Gx1[73]*Gx2[78] + Gx1[87]*Gx2[92] + Gx1[101]*Gx2[106] + Gx1[115]*Gx2[120] + Gx1[129]*Gx2[134] + Gx1[143]*Gx2[148] + Gx1[157]*Gx2[162] + Gx1[171]*Gx2[176] + Gx1[185]*Gx2[190];
nmpcWorkspace.H[651] += + Gx1[3]*Gx2[9] + Gx1[17]*Gx2[23] + Gx1[31]*Gx2[37] + Gx1[45]*Gx2[51] + Gx1[59]*Gx2[65] + Gx1[73]*Gx2[79] + Gx1[87]*Gx2[93] + Gx1[101]*Gx2[107] + Gx1[115]*Gx2[121] + Gx1[129]*Gx2[135] + Gx1[143]*Gx2[149] + Gx1[157]*Gx2[163] + Gx1[171]*Gx2[177] + Gx1[185]*Gx2[191];
nmpcWorkspace.H[652] += + Gx1[3]*Gx2[10] + Gx1[17]*Gx2[24] + Gx1[31]*Gx2[38] + Gx1[45]*Gx2[52] + Gx1[59]*Gx2[66] + Gx1[73]*Gx2[80] + Gx1[87]*Gx2[94] + Gx1[101]*Gx2[108] + Gx1[115]*Gx2[122] + Gx1[129]*Gx2[136] + Gx1[143]*Gx2[150] + Gx1[157]*Gx2[164] + Gx1[171]*Gx2[178] + Gx1[185]*Gx2[192];
nmpcWorkspace.H[653] += + Gx1[3]*Gx2[11] + Gx1[17]*Gx2[25] + Gx1[31]*Gx2[39] + Gx1[45]*Gx2[53] + Gx1[59]*Gx2[67] + Gx1[73]*Gx2[81] + Gx1[87]*Gx2[95] + Gx1[101]*Gx2[109] + Gx1[115]*Gx2[123] + Gx1[129]*Gx2[137] + Gx1[143]*Gx2[151] + Gx1[157]*Gx2[165] + Gx1[171]*Gx2[179] + Gx1[185]*Gx2[193];
nmpcWorkspace.H[654] += + Gx1[3]*Gx2[12] + Gx1[17]*Gx2[26] + Gx1[31]*Gx2[40] + Gx1[45]*Gx2[54] + Gx1[59]*Gx2[68] + Gx1[73]*Gx2[82] + Gx1[87]*Gx2[96] + Gx1[101]*Gx2[110] + Gx1[115]*Gx2[124] + Gx1[129]*Gx2[138] + Gx1[143]*Gx2[152] + Gx1[157]*Gx2[166] + Gx1[171]*Gx2[180] + Gx1[185]*Gx2[194];
nmpcWorkspace.H[655] += + Gx1[3]*Gx2[13] + Gx1[17]*Gx2[27] + Gx1[31]*Gx2[41] + Gx1[45]*Gx2[55] + Gx1[59]*Gx2[69] + Gx1[73]*Gx2[83] + Gx1[87]*Gx2[97] + Gx1[101]*Gx2[111] + Gx1[115]*Gx2[125] + Gx1[129]*Gx2[139] + Gx1[143]*Gx2[153] + Gx1[157]*Gx2[167] + Gx1[171]*Gx2[181] + Gx1[185]*Gx2[195];
nmpcWorkspace.H[856] += + Gx1[4]*Gx2[0] + Gx1[18]*Gx2[14] + Gx1[32]*Gx2[28] + Gx1[46]*Gx2[42] + Gx1[60]*Gx2[56] + Gx1[74]*Gx2[70] + Gx1[88]*Gx2[84] + Gx1[102]*Gx2[98] + Gx1[116]*Gx2[112] + Gx1[130]*Gx2[126] + Gx1[144]*Gx2[140] + Gx1[158]*Gx2[154] + Gx1[172]*Gx2[168] + Gx1[186]*Gx2[182];
nmpcWorkspace.H[857] += + Gx1[4]*Gx2[1] + Gx1[18]*Gx2[15] + Gx1[32]*Gx2[29] + Gx1[46]*Gx2[43] + Gx1[60]*Gx2[57] + Gx1[74]*Gx2[71] + Gx1[88]*Gx2[85] + Gx1[102]*Gx2[99] + Gx1[116]*Gx2[113] + Gx1[130]*Gx2[127] + Gx1[144]*Gx2[141] + Gx1[158]*Gx2[155] + Gx1[172]*Gx2[169] + Gx1[186]*Gx2[183];
nmpcWorkspace.H[858] += + Gx1[4]*Gx2[2] + Gx1[18]*Gx2[16] + Gx1[32]*Gx2[30] + Gx1[46]*Gx2[44] + Gx1[60]*Gx2[58] + Gx1[74]*Gx2[72] + Gx1[88]*Gx2[86] + Gx1[102]*Gx2[100] + Gx1[116]*Gx2[114] + Gx1[130]*Gx2[128] + Gx1[144]*Gx2[142] + Gx1[158]*Gx2[156] + Gx1[172]*Gx2[170] + Gx1[186]*Gx2[184];
nmpcWorkspace.H[859] += + Gx1[4]*Gx2[3] + Gx1[18]*Gx2[17] + Gx1[32]*Gx2[31] + Gx1[46]*Gx2[45] + Gx1[60]*Gx2[59] + Gx1[74]*Gx2[73] + Gx1[88]*Gx2[87] + Gx1[102]*Gx2[101] + Gx1[116]*Gx2[115] + Gx1[130]*Gx2[129] + Gx1[144]*Gx2[143] + Gx1[158]*Gx2[157] + Gx1[172]*Gx2[171] + Gx1[186]*Gx2[185];
nmpcWorkspace.H[860] += + Gx1[4]*Gx2[4] + Gx1[18]*Gx2[18] + Gx1[32]*Gx2[32] + Gx1[46]*Gx2[46] + Gx1[60]*Gx2[60] + Gx1[74]*Gx2[74] + Gx1[88]*Gx2[88] + Gx1[102]*Gx2[102] + Gx1[116]*Gx2[116] + Gx1[130]*Gx2[130] + Gx1[144]*Gx2[144] + Gx1[158]*Gx2[158] + Gx1[172]*Gx2[172] + Gx1[186]*Gx2[186];
nmpcWorkspace.H[861] += + Gx1[4]*Gx2[5] + Gx1[18]*Gx2[19] + Gx1[32]*Gx2[33] + Gx1[46]*Gx2[47] + Gx1[60]*Gx2[61] + Gx1[74]*Gx2[75] + Gx1[88]*Gx2[89] + Gx1[102]*Gx2[103] + Gx1[116]*Gx2[117] + Gx1[130]*Gx2[131] + Gx1[144]*Gx2[145] + Gx1[158]*Gx2[159] + Gx1[172]*Gx2[173] + Gx1[186]*Gx2[187];
nmpcWorkspace.H[862] += + Gx1[4]*Gx2[6] + Gx1[18]*Gx2[20] + Gx1[32]*Gx2[34] + Gx1[46]*Gx2[48] + Gx1[60]*Gx2[62] + Gx1[74]*Gx2[76] + Gx1[88]*Gx2[90] + Gx1[102]*Gx2[104] + Gx1[116]*Gx2[118] + Gx1[130]*Gx2[132] + Gx1[144]*Gx2[146] + Gx1[158]*Gx2[160] + Gx1[172]*Gx2[174] + Gx1[186]*Gx2[188];
nmpcWorkspace.H[863] += + Gx1[4]*Gx2[7] + Gx1[18]*Gx2[21] + Gx1[32]*Gx2[35] + Gx1[46]*Gx2[49] + Gx1[60]*Gx2[63] + Gx1[74]*Gx2[77] + Gx1[88]*Gx2[91] + Gx1[102]*Gx2[105] + Gx1[116]*Gx2[119] + Gx1[130]*Gx2[133] + Gx1[144]*Gx2[147] + Gx1[158]*Gx2[161] + Gx1[172]*Gx2[175] + Gx1[186]*Gx2[189];
nmpcWorkspace.H[864] += + Gx1[4]*Gx2[8] + Gx1[18]*Gx2[22] + Gx1[32]*Gx2[36] + Gx1[46]*Gx2[50] + Gx1[60]*Gx2[64] + Gx1[74]*Gx2[78] + Gx1[88]*Gx2[92] + Gx1[102]*Gx2[106] + Gx1[116]*Gx2[120] + Gx1[130]*Gx2[134] + Gx1[144]*Gx2[148] + Gx1[158]*Gx2[162] + Gx1[172]*Gx2[176] + Gx1[186]*Gx2[190];
nmpcWorkspace.H[865] += + Gx1[4]*Gx2[9] + Gx1[18]*Gx2[23] + Gx1[32]*Gx2[37] + Gx1[46]*Gx2[51] + Gx1[60]*Gx2[65] + Gx1[74]*Gx2[79] + Gx1[88]*Gx2[93] + Gx1[102]*Gx2[107] + Gx1[116]*Gx2[121] + Gx1[130]*Gx2[135] + Gx1[144]*Gx2[149] + Gx1[158]*Gx2[163] + Gx1[172]*Gx2[177] + Gx1[186]*Gx2[191];
nmpcWorkspace.H[866] += + Gx1[4]*Gx2[10] + Gx1[18]*Gx2[24] + Gx1[32]*Gx2[38] + Gx1[46]*Gx2[52] + Gx1[60]*Gx2[66] + Gx1[74]*Gx2[80] + Gx1[88]*Gx2[94] + Gx1[102]*Gx2[108] + Gx1[116]*Gx2[122] + Gx1[130]*Gx2[136] + Gx1[144]*Gx2[150] + Gx1[158]*Gx2[164] + Gx1[172]*Gx2[178] + Gx1[186]*Gx2[192];
nmpcWorkspace.H[867] += + Gx1[4]*Gx2[11] + Gx1[18]*Gx2[25] + Gx1[32]*Gx2[39] + Gx1[46]*Gx2[53] + Gx1[60]*Gx2[67] + Gx1[74]*Gx2[81] + Gx1[88]*Gx2[95] + Gx1[102]*Gx2[109] + Gx1[116]*Gx2[123] + Gx1[130]*Gx2[137] + Gx1[144]*Gx2[151] + Gx1[158]*Gx2[165] + Gx1[172]*Gx2[179] + Gx1[186]*Gx2[193];
nmpcWorkspace.H[868] += + Gx1[4]*Gx2[12] + Gx1[18]*Gx2[26] + Gx1[32]*Gx2[40] + Gx1[46]*Gx2[54] + Gx1[60]*Gx2[68] + Gx1[74]*Gx2[82] + Gx1[88]*Gx2[96] + Gx1[102]*Gx2[110] + Gx1[116]*Gx2[124] + Gx1[130]*Gx2[138] + Gx1[144]*Gx2[152] + Gx1[158]*Gx2[166] + Gx1[172]*Gx2[180] + Gx1[186]*Gx2[194];
nmpcWorkspace.H[869] += + Gx1[4]*Gx2[13] + Gx1[18]*Gx2[27] + Gx1[32]*Gx2[41] + Gx1[46]*Gx2[55] + Gx1[60]*Gx2[69] + Gx1[74]*Gx2[83] + Gx1[88]*Gx2[97] + Gx1[102]*Gx2[111] + Gx1[116]*Gx2[125] + Gx1[130]*Gx2[139] + Gx1[144]*Gx2[153] + Gx1[158]*Gx2[167] + Gx1[172]*Gx2[181] + Gx1[186]*Gx2[195];
nmpcWorkspace.H[1070] += + Gx1[5]*Gx2[0] + Gx1[19]*Gx2[14] + Gx1[33]*Gx2[28] + Gx1[47]*Gx2[42] + Gx1[61]*Gx2[56] + Gx1[75]*Gx2[70] + Gx1[89]*Gx2[84] + Gx1[103]*Gx2[98] + Gx1[117]*Gx2[112] + Gx1[131]*Gx2[126] + Gx1[145]*Gx2[140] + Gx1[159]*Gx2[154] + Gx1[173]*Gx2[168] + Gx1[187]*Gx2[182];
nmpcWorkspace.H[1071] += + Gx1[5]*Gx2[1] + Gx1[19]*Gx2[15] + Gx1[33]*Gx2[29] + Gx1[47]*Gx2[43] + Gx1[61]*Gx2[57] + Gx1[75]*Gx2[71] + Gx1[89]*Gx2[85] + Gx1[103]*Gx2[99] + Gx1[117]*Gx2[113] + Gx1[131]*Gx2[127] + Gx1[145]*Gx2[141] + Gx1[159]*Gx2[155] + Gx1[173]*Gx2[169] + Gx1[187]*Gx2[183];
nmpcWorkspace.H[1072] += + Gx1[5]*Gx2[2] + Gx1[19]*Gx2[16] + Gx1[33]*Gx2[30] + Gx1[47]*Gx2[44] + Gx1[61]*Gx2[58] + Gx1[75]*Gx2[72] + Gx1[89]*Gx2[86] + Gx1[103]*Gx2[100] + Gx1[117]*Gx2[114] + Gx1[131]*Gx2[128] + Gx1[145]*Gx2[142] + Gx1[159]*Gx2[156] + Gx1[173]*Gx2[170] + Gx1[187]*Gx2[184];
nmpcWorkspace.H[1073] += + Gx1[5]*Gx2[3] + Gx1[19]*Gx2[17] + Gx1[33]*Gx2[31] + Gx1[47]*Gx2[45] + Gx1[61]*Gx2[59] + Gx1[75]*Gx2[73] + Gx1[89]*Gx2[87] + Gx1[103]*Gx2[101] + Gx1[117]*Gx2[115] + Gx1[131]*Gx2[129] + Gx1[145]*Gx2[143] + Gx1[159]*Gx2[157] + Gx1[173]*Gx2[171] + Gx1[187]*Gx2[185];
nmpcWorkspace.H[1074] += + Gx1[5]*Gx2[4] + Gx1[19]*Gx2[18] + Gx1[33]*Gx2[32] + Gx1[47]*Gx2[46] + Gx1[61]*Gx2[60] + Gx1[75]*Gx2[74] + Gx1[89]*Gx2[88] + Gx1[103]*Gx2[102] + Gx1[117]*Gx2[116] + Gx1[131]*Gx2[130] + Gx1[145]*Gx2[144] + Gx1[159]*Gx2[158] + Gx1[173]*Gx2[172] + Gx1[187]*Gx2[186];
nmpcWorkspace.H[1075] += + Gx1[5]*Gx2[5] + Gx1[19]*Gx2[19] + Gx1[33]*Gx2[33] + Gx1[47]*Gx2[47] + Gx1[61]*Gx2[61] + Gx1[75]*Gx2[75] + Gx1[89]*Gx2[89] + Gx1[103]*Gx2[103] + Gx1[117]*Gx2[117] + Gx1[131]*Gx2[131] + Gx1[145]*Gx2[145] + Gx1[159]*Gx2[159] + Gx1[173]*Gx2[173] + Gx1[187]*Gx2[187];
nmpcWorkspace.H[1076] += + Gx1[5]*Gx2[6] + Gx1[19]*Gx2[20] + Gx1[33]*Gx2[34] + Gx1[47]*Gx2[48] + Gx1[61]*Gx2[62] + Gx1[75]*Gx2[76] + Gx1[89]*Gx2[90] + Gx1[103]*Gx2[104] + Gx1[117]*Gx2[118] + Gx1[131]*Gx2[132] + Gx1[145]*Gx2[146] + Gx1[159]*Gx2[160] + Gx1[173]*Gx2[174] + Gx1[187]*Gx2[188];
nmpcWorkspace.H[1077] += + Gx1[5]*Gx2[7] + Gx1[19]*Gx2[21] + Gx1[33]*Gx2[35] + Gx1[47]*Gx2[49] + Gx1[61]*Gx2[63] + Gx1[75]*Gx2[77] + Gx1[89]*Gx2[91] + Gx1[103]*Gx2[105] + Gx1[117]*Gx2[119] + Gx1[131]*Gx2[133] + Gx1[145]*Gx2[147] + Gx1[159]*Gx2[161] + Gx1[173]*Gx2[175] + Gx1[187]*Gx2[189];
nmpcWorkspace.H[1078] += + Gx1[5]*Gx2[8] + Gx1[19]*Gx2[22] + Gx1[33]*Gx2[36] + Gx1[47]*Gx2[50] + Gx1[61]*Gx2[64] + Gx1[75]*Gx2[78] + Gx1[89]*Gx2[92] + Gx1[103]*Gx2[106] + Gx1[117]*Gx2[120] + Gx1[131]*Gx2[134] + Gx1[145]*Gx2[148] + Gx1[159]*Gx2[162] + Gx1[173]*Gx2[176] + Gx1[187]*Gx2[190];
nmpcWorkspace.H[1079] += + Gx1[5]*Gx2[9] + Gx1[19]*Gx2[23] + Gx1[33]*Gx2[37] + Gx1[47]*Gx2[51] + Gx1[61]*Gx2[65] + Gx1[75]*Gx2[79] + Gx1[89]*Gx2[93] + Gx1[103]*Gx2[107] + Gx1[117]*Gx2[121] + Gx1[131]*Gx2[135] + Gx1[145]*Gx2[149] + Gx1[159]*Gx2[163] + Gx1[173]*Gx2[177] + Gx1[187]*Gx2[191];
nmpcWorkspace.H[1080] += + Gx1[5]*Gx2[10] + Gx1[19]*Gx2[24] + Gx1[33]*Gx2[38] + Gx1[47]*Gx2[52] + Gx1[61]*Gx2[66] + Gx1[75]*Gx2[80] + Gx1[89]*Gx2[94] + Gx1[103]*Gx2[108] + Gx1[117]*Gx2[122] + Gx1[131]*Gx2[136] + Gx1[145]*Gx2[150] + Gx1[159]*Gx2[164] + Gx1[173]*Gx2[178] + Gx1[187]*Gx2[192];
nmpcWorkspace.H[1081] += + Gx1[5]*Gx2[11] + Gx1[19]*Gx2[25] + Gx1[33]*Gx2[39] + Gx1[47]*Gx2[53] + Gx1[61]*Gx2[67] + Gx1[75]*Gx2[81] + Gx1[89]*Gx2[95] + Gx1[103]*Gx2[109] + Gx1[117]*Gx2[123] + Gx1[131]*Gx2[137] + Gx1[145]*Gx2[151] + Gx1[159]*Gx2[165] + Gx1[173]*Gx2[179] + Gx1[187]*Gx2[193];
nmpcWorkspace.H[1082] += + Gx1[5]*Gx2[12] + Gx1[19]*Gx2[26] + Gx1[33]*Gx2[40] + Gx1[47]*Gx2[54] + Gx1[61]*Gx2[68] + Gx1[75]*Gx2[82] + Gx1[89]*Gx2[96] + Gx1[103]*Gx2[110] + Gx1[117]*Gx2[124] + Gx1[131]*Gx2[138] + Gx1[145]*Gx2[152] + Gx1[159]*Gx2[166] + Gx1[173]*Gx2[180] + Gx1[187]*Gx2[194];
nmpcWorkspace.H[1083] += + Gx1[5]*Gx2[13] + Gx1[19]*Gx2[27] + Gx1[33]*Gx2[41] + Gx1[47]*Gx2[55] + Gx1[61]*Gx2[69] + Gx1[75]*Gx2[83] + Gx1[89]*Gx2[97] + Gx1[103]*Gx2[111] + Gx1[117]*Gx2[125] + Gx1[131]*Gx2[139] + Gx1[145]*Gx2[153] + Gx1[159]*Gx2[167] + Gx1[173]*Gx2[181] + Gx1[187]*Gx2[195];
nmpcWorkspace.H[1284] += + Gx1[6]*Gx2[0] + Gx1[20]*Gx2[14] + Gx1[34]*Gx2[28] + Gx1[48]*Gx2[42] + Gx1[62]*Gx2[56] + Gx1[76]*Gx2[70] + Gx1[90]*Gx2[84] + Gx1[104]*Gx2[98] + Gx1[118]*Gx2[112] + Gx1[132]*Gx2[126] + Gx1[146]*Gx2[140] + Gx1[160]*Gx2[154] + Gx1[174]*Gx2[168] + Gx1[188]*Gx2[182];
nmpcWorkspace.H[1285] += + Gx1[6]*Gx2[1] + Gx1[20]*Gx2[15] + Gx1[34]*Gx2[29] + Gx1[48]*Gx2[43] + Gx1[62]*Gx2[57] + Gx1[76]*Gx2[71] + Gx1[90]*Gx2[85] + Gx1[104]*Gx2[99] + Gx1[118]*Gx2[113] + Gx1[132]*Gx2[127] + Gx1[146]*Gx2[141] + Gx1[160]*Gx2[155] + Gx1[174]*Gx2[169] + Gx1[188]*Gx2[183];
nmpcWorkspace.H[1286] += + Gx1[6]*Gx2[2] + Gx1[20]*Gx2[16] + Gx1[34]*Gx2[30] + Gx1[48]*Gx2[44] + Gx1[62]*Gx2[58] + Gx1[76]*Gx2[72] + Gx1[90]*Gx2[86] + Gx1[104]*Gx2[100] + Gx1[118]*Gx2[114] + Gx1[132]*Gx2[128] + Gx1[146]*Gx2[142] + Gx1[160]*Gx2[156] + Gx1[174]*Gx2[170] + Gx1[188]*Gx2[184];
nmpcWorkspace.H[1287] += + Gx1[6]*Gx2[3] + Gx1[20]*Gx2[17] + Gx1[34]*Gx2[31] + Gx1[48]*Gx2[45] + Gx1[62]*Gx2[59] + Gx1[76]*Gx2[73] + Gx1[90]*Gx2[87] + Gx1[104]*Gx2[101] + Gx1[118]*Gx2[115] + Gx1[132]*Gx2[129] + Gx1[146]*Gx2[143] + Gx1[160]*Gx2[157] + Gx1[174]*Gx2[171] + Gx1[188]*Gx2[185];
nmpcWorkspace.H[1288] += + Gx1[6]*Gx2[4] + Gx1[20]*Gx2[18] + Gx1[34]*Gx2[32] + Gx1[48]*Gx2[46] + Gx1[62]*Gx2[60] + Gx1[76]*Gx2[74] + Gx1[90]*Gx2[88] + Gx1[104]*Gx2[102] + Gx1[118]*Gx2[116] + Gx1[132]*Gx2[130] + Gx1[146]*Gx2[144] + Gx1[160]*Gx2[158] + Gx1[174]*Gx2[172] + Gx1[188]*Gx2[186];
nmpcWorkspace.H[1289] += + Gx1[6]*Gx2[5] + Gx1[20]*Gx2[19] + Gx1[34]*Gx2[33] + Gx1[48]*Gx2[47] + Gx1[62]*Gx2[61] + Gx1[76]*Gx2[75] + Gx1[90]*Gx2[89] + Gx1[104]*Gx2[103] + Gx1[118]*Gx2[117] + Gx1[132]*Gx2[131] + Gx1[146]*Gx2[145] + Gx1[160]*Gx2[159] + Gx1[174]*Gx2[173] + Gx1[188]*Gx2[187];
nmpcWorkspace.H[1290] += + Gx1[6]*Gx2[6] + Gx1[20]*Gx2[20] + Gx1[34]*Gx2[34] + Gx1[48]*Gx2[48] + Gx1[62]*Gx2[62] + Gx1[76]*Gx2[76] + Gx1[90]*Gx2[90] + Gx1[104]*Gx2[104] + Gx1[118]*Gx2[118] + Gx1[132]*Gx2[132] + Gx1[146]*Gx2[146] + Gx1[160]*Gx2[160] + Gx1[174]*Gx2[174] + Gx1[188]*Gx2[188];
nmpcWorkspace.H[1291] += + Gx1[6]*Gx2[7] + Gx1[20]*Gx2[21] + Gx1[34]*Gx2[35] + Gx1[48]*Gx2[49] + Gx1[62]*Gx2[63] + Gx1[76]*Gx2[77] + Gx1[90]*Gx2[91] + Gx1[104]*Gx2[105] + Gx1[118]*Gx2[119] + Gx1[132]*Gx2[133] + Gx1[146]*Gx2[147] + Gx1[160]*Gx2[161] + Gx1[174]*Gx2[175] + Gx1[188]*Gx2[189];
nmpcWorkspace.H[1292] += + Gx1[6]*Gx2[8] + Gx1[20]*Gx2[22] + Gx1[34]*Gx2[36] + Gx1[48]*Gx2[50] + Gx1[62]*Gx2[64] + Gx1[76]*Gx2[78] + Gx1[90]*Gx2[92] + Gx1[104]*Gx2[106] + Gx1[118]*Gx2[120] + Gx1[132]*Gx2[134] + Gx1[146]*Gx2[148] + Gx1[160]*Gx2[162] + Gx1[174]*Gx2[176] + Gx1[188]*Gx2[190];
nmpcWorkspace.H[1293] += + Gx1[6]*Gx2[9] + Gx1[20]*Gx2[23] + Gx1[34]*Gx2[37] + Gx1[48]*Gx2[51] + Gx1[62]*Gx2[65] + Gx1[76]*Gx2[79] + Gx1[90]*Gx2[93] + Gx1[104]*Gx2[107] + Gx1[118]*Gx2[121] + Gx1[132]*Gx2[135] + Gx1[146]*Gx2[149] + Gx1[160]*Gx2[163] + Gx1[174]*Gx2[177] + Gx1[188]*Gx2[191];
nmpcWorkspace.H[1294] += + Gx1[6]*Gx2[10] + Gx1[20]*Gx2[24] + Gx1[34]*Gx2[38] + Gx1[48]*Gx2[52] + Gx1[62]*Gx2[66] + Gx1[76]*Gx2[80] + Gx1[90]*Gx2[94] + Gx1[104]*Gx2[108] + Gx1[118]*Gx2[122] + Gx1[132]*Gx2[136] + Gx1[146]*Gx2[150] + Gx1[160]*Gx2[164] + Gx1[174]*Gx2[178] + Gx1[188]*Gx2[192];
nmpcWorkspace.H[1295] += + Gx1[6]*Gx2[11] + Gx1[20]*Gx2[25] + Gx1[34]*Gx2[39] + Gx1[48]*Gx2[53] + Gx1[62]*Gx2[67] + Gx1[76]*Gx2[81] + Gx1[90]*Gx2[95] + Gx1[104]*Gx2[109] + Gx1[118]*Gx2[123] + Gx1[132]*Gx2[137] + Gx1[146]*Gx2[151] + Gx1[160]*Gx2[165] + Gx1[174]*Gx2[179] + Gx1[188]*Gx2[193];
nmpcWorkspace.H[1296] += + Gx1[6]*Gx2[12] + Gx1[20]*Gx2[26] + Gx1[34]*Gx2[40] + Gx1[48]*Gx2[54] + Gx1[62]*Gx2[68] + Gx1[76]*Gx2[82] + Gx1[90]*Gx2[96] + Gx1[104]*Gx2[110] + Gx1[118]*Gx2[124] + Gx1[132]*Gx2[138] + Gx1[146]*Gx2[152] + Gx1[160]*Gx2[166] + Gx1[174]*Gx2[180] + Gx1[188]*Gx2[194];
nmpcWorkspace.H[1297] += + Gx1[6]*Gx2[13] + Gx1[20]*Gx2[27] + Gx1[34]*Gx2[41] + Gx1[48]*Gx2[55] + Gx1[62]*Gx2[69] + Gx1[76]*Gx2[83] + Gx1[90]*Gx2[97] + Gx1[104]*Gx2[111] + Gx1[118]*Gx2[125] + Gx1[132]*Gx2[139] + Gx1[146]*Gx2[153] + Gx1[160]*Gx2[167] + Gx1[174]*Gx2[181] + Gx1[188]*Gx2[195];
nmpcWorkspace.H[1498] += + Gx1[7]*Gx2[0] + Gx1[21]*Gx2[14] + Gx1[35]*Gx2[28] + Gx1[49]*Gx2[42] + Gx1[63]*Gx2[56] + Gx1[77]*Gx2[70] + Gx1[91]*Gx2[84] + Gx1[105]*Gx2[98] + Gx1[119]*Gx2[112] + Gx1[133]*Gx2[126] + Gx1[147]*Gx2[140] + Gx1[161]*Gx2[154] + Gx1[175]*Gx2[168] + Gx1[189]*Gx2[182];
nmpcWorkspace.H[1499] += + Gx1[7]*Gx2[1] + Gx1[21]*Gx2[15] + Gx1[35]*Gx2[29] + Gx1[49]*Gx2[43] + Gx1[63]*Gx2[57] + Gx1[77]*Gx2[71] + Gx1[91]*Gx2[85] + Gx1[105]*Gx2[99] + Gx1[119]*Gx2[113] + Gx1[133]*Gx2[127] + Gx1[147]*Gx2[141] + Gx1[161]*Gx2[155] + Gx1[175]*Gx2[169] + Gx1[189]*Gx2[183];
nmpcWorkspace.H[1500] += + Gx1[7]*Gx2[2] + Gx1[21]*Gx2[16] + Gx1[35]*Gx2[30] + Gx1[49]*Gx2[44] + Gx1[63]*Gx2[58] + Gx1[77]*Gx2[72] + Gx1[91]*Gx2[86] + Gx1[105]*Gx2[100] + Gx1[119]*Gx2[114] + Gx1[133]*Gx2[128] + Gx1[147]*Gx2[142] + Gx1[161]*Gx2[156] + Gx1[175]*Gx2[170] + Gx1[189]*Gx2[184];
nmpcWorkspace.H[1501] += + Gx1[7]*Gx2[3] + Gx1[21]*Gx2[17] + Gx1[35]*Gx2[31] + Gx1[49]*Gx2[45] + Gx1[63]*Gx2[59] + Gx1[77]*Gx2[73] + Gx1[91]*Gx2[87] + Gx1[105]*Gx2[101] + Gx1[119]*Gx2[115] + Gx1[133]*Gx2[129] + Gx1[147]*Gx2[143] + Gx1[161]*Gx2[157] + Gx1[175]*Gx2[171] + Gx1[189]*Gx2[185];
nmpcWorkspace.H[1502] += + Gx1[7]*Gx2[4] + Gx1[21]*Gx2[18] + Gx1[35]*Gx2[32] + Gx1[49]*Gx2[46] + Gx1[63]*Gx2[60] + Gx1[77]*Gx2[74] + Gx1[91]*Gx2[88] + Gx1[105]*Gx2[102] + Gx1[119]*Gx2[116] + Gx1[133]*Gx2[130] + Gx1[147]*Gx2[144] + Gx1[161]*Gx2[158] + Gx1[175]*Gx2[172] + Gx1[189]*Gx2[186];
nmpcWorkspace.H[1503] += + Gx1[7]*Gx2[5] + Gx1[21]*Gx2[19] + Gx1[35]*Gx2[33] + Gx1[49]*Gx2[47] + Gx1[63]*Gx2[61] + Gx1[77]*Gx2[75] + Gx1[91]*Gx2[89] + Gx1[105]*Gx2[103] + Gx1[119]*Gx2[117] + Gx1[133]*Gx2[131] + Gx1[147]*Gx2[145] + Gx1[161]*Gx2[159] + Gx1[175]*Gx2[173] + Gx1[189]*Gx2[187];
nmpcWorkspace.H[1504] += + Gx1[7]*Gx2[6] + Gx1[21]*Gx2[20] + Gx1[35]*Gx2[34] + Gx1[49]*Gx2[48] + Gx1[63]*Gx2[62] + Gx1[77]*Gx2[76] + Gx1[91]*Gx2[90] + Gx1[105]*Gx2[104] + Gx1[119]*Gx2[118] + Gx1[133]*Gx2[132] + Gx1[147]*Gx2[146] + Gx1[161]*Gx2[160] + Gx1[175]*Gx2[174] + Gx1[189]*Gx2[188];
nmpcWorkspace.H[1505] += + Gx1[7]*Gx2[7] + Gx1[21]*Gx2[21] + Gx1[35]*Gx2[35] + Gx1[49]*Gx2[49] + Gx1[63]*Gx2[63] + Gx1[77]*Gx2[77] + Gx1[91]*Gx2[91] + Gx1[105]*Gx2[105] + Gx1[119]*Gx2[119] + Gx1[133]*Gx2[133] + Gx1[147]*Gx2[147] + Gx1[161]*Gx2[161] + Gx1[175]*Gx2[175] + Gx1[189]*Gx2[189];
nmpcWorkspace.H[1506] += + Gx1[7]*Gx2[8] + Gx1[21]*Gx2[22] + Gx1[35]*Gx2[36] + Gx1[49]*Gx2[50] + Gx1[63]*Gx2[64] + Gx1[77]*Gx2[78] + Gx1[91]*Gx2[92] + Gx1[105]*Gx2[106] + Gx1[119]*Gx2[120] + Gx1[133]*Gx2[134] + Gx1[147]*Gx2[148] + Gx1[161]*Gx2[162] + Gx1[175]*Gx2[176] + Gx1[189]*Gx2[190];
nmpcWorkspace.H[1507] += + Gx1[7]*Gx2[9] + Gx1[21]*Gx2[23] + Gx1[35]*Gx2[37] + Gx1[49]*Gx2[51] + Gx1[63]*Gx2[65] + Gx1[77]*Gx2[79] + Gx1[91]*Gx2[93] + Gx1[105]*Gx2[107] + Gx1[119]*Gx2[121] + Gx1[133]*Gx2[135] + Gx1[147]*Gx2[149] + Gx1[161]*Gx2[163] + Gx1[175]*Gx2[177] + Gx1[189]*Gx2[191];
nmpcWorkspace.H[1508] += + Gx1[7]*Gx2[10] + Gx1[21]*Gx2[24] + Gx1[35]*Gx2[38] + Gx1[49]*Gx2[52] + Gx1[63]*Gx2[66] + Gx1[77]*Gx2[80] + Gx1[91]*Gx2[94] + Gx1[105]*Gx2[108] + Gx1[119]*Gx2[122] + Gx1[133]*Gx2[136] + Gx1[147]*Gx2[150] + Gx1[161]*Gx2[164] + Gx1[175]*Gx2[178] + Gx1[189]*Gx2[192];
nmpcWorkspace.H[1509] += + Gx1[7]*Gx2[11] + Gx1[21]*Gx2[25] + Gx1[35]*Gx2[39] + Gx1[49]*Gx2[53] + Gx1[63]*Gx2[67] + Gx1[77]*Gx2[81] + Gx1[91]*Gx2[95] + Gx1[105]*Gx2[109] + Gx1[119]*Gx2[123] + Gx1[133]*Gx2[137] + Gx1[147]*Gx2[151] + Gx1[161]*Gx2[165] + Gx1[175]*Gx2[179] + Gx1[189]*Gx2[193];
nmpcWorkspace.H[1510] += + Gx1[7]*Gx2[12] + Gx1[21]*Gx2[26] + Gx1[35]*Gx2[40] + Gx1[49]*Gx2[54] + Gx1[63]*Gx2[68] + Gx1[77]*Gx2[82] + Gx1[91]*Gx2[96] + Gx1[105]*Gx2[110] + Gx1[119]*Gx2[124] + Gx1[133]*Gx2[138] + Gx1[147]*Gx2[152] + Gx1[161]*Gx2[166] + Gx1[175]*Gx2[180] + Gx1[189]*Gx2[194];
nmpcWorkspace.H[1511] += + Gx1[7]*Gx2[13] + Gx1[21]*Gx2[27] + Gx1[35]*Gx2[41] + Gx1[49]*Gx2[55] + Gx1[63]*Gx2[69] + Gx1[77]*Gx2[83] + Gx1[91]*Gx2[97] + Gx1[105]*Gx2[111] + Gx1[119]*Gx2[125] + Gx1[133]*Gx2[139] + Gx1[147]*Gx2[153] + Gx1[161]*Gx2[167] + Gx1[175]*Gx2[181] + Gx1[189]*Gx2[195];
nmpcWorkspace.H[1712] += + Gx1[8]*Gx2[0] + Gx1[22]*Gx2[14] + Gx1[36]*Gx2[28] + Gx1[50]*Gx2[42] + Gx1[64]*Gx2[56] + Gx1[78]*Gx2[70] + Gx1[92]*Gx2[84] + Gx1[106]*Gx2[98] + Gx1[120]*Gx2[112] + Gx1[134]*Gx2[126] + Gx1[148]*Gx2[140] + Gx1[162]*Gx2[154] + Gx1[176]*Gx2[168] + Gx1[190]*Gx2[182];
nmpcWorkspace.H[1713] += + Gx1[8]*Gx2[1] + Gx1[22]*Gx2[15] + Gx1[36]*Gx2[29] + Gx1[50]*Gx2[43] + Gx1[64]*Gx2[57] + Gx1[78]*Gx2[71] + Gx1[92]*Gx2[85] + Gx1[106]*Gx2[99] + Gx1[120]*Gx2[113] + Gx1[134]*Gx2[127] + Gx1[148]*Gx2[141] + Gx1[162]*Gx2[155] + Gx1[176]*Gx2[169] + Gx1[190]*Gx2[183];
nmpcWorkspace.H[1714] += + Gx1[8]*Gx2[2] + Gx1[22]*Gx2[16] + Gx1[36]*Gx2[30] + Gx1[50]*Gx2[44] + Gx1[64]*Gx2[58] + Gx1[78]*Gx2[72] + Gx1[92]*Gx2[86] + Gx1[106]*Gx2[100] + Gx1[120]*Gx2[114] + Gx1[134]*Gx2[128] + Gx1[148]*Gx2[142] + Gx1[162]*Gx2[156] + Gx1[176]*Gx2[170] + Gx1[190]*Gx2[184];
nmpcWorkspace.H[1715] += + Gx1[8]*Gx2[3] + Gx1[22]*Gx2[17] + Gx1[36]*Gx2[31] + Gx1[50]*Gx2[45] + Gx1[64]*Gx2[59] + Gx1[78]*Gx2[73] + Gx1[92]*Gx2[87] + Gx1[106]*Gx2[101] + Gx1[120]*Gx2[115] + Gx1[134]*Gx2[129] + Gx1[148]*Gx2[143] + Gx1[162]*Gx2[157] + Gx1[176]*Gx2[171] + Gx1[190]*Gx2[185];
nmpcWorkspace.H[1716] += + Gx1[8]*Gx2[4] + Gx1[22]*Gx2[18] + Gx1[36]*Gx2[32] + Gx1[50]*Gx2[46] + Gx1[64]*Gx2[60] + Gx1[78]*Gx2[74] + Gx1[92]*Gx2[88] + Gx1[106]*Gx2[102] + Gx1[120]*Gx2[116] + Gx1[134]*Gx2[130] + Gx1[148]*Gx2[144] + Gx1[162]*Gx2[158] + Gx1[176]*Gx2[172] + Gx1[190]*Gx2[186];
nmpcWorkspace.H[1717] += + Gx1[8]*Gx2[5] + Gx1[22]*Gx2[19] + Gx1[36]*Gx2[33] + Gx1[50]*Gx2[47] + Gx1[64]*Gx2[61] + Gx1[78]*Gx2[75] + Gx1[92]*Gx2[89] + Gx1[106]*Gx2[103] + Gx1[120]*Gx2[117] + Gx1[134]*Gx2[131] + Gx1[148]*Gx2[145] + Gx1[162]*Gx2[159] + Gx1[176]*Gx2[173] + Gx1[190]*Gx2[187];
nmpcWorkspace.H[1718] += + Gx1[8]*Gx2[6] + Gx1[22]*Gx2[20] + Gx1[36]*Gx2[34] + Gx1[50]*Gx2[48] + Gx1[64]*Gx2[62] + Gx1[78]*Gx2[76] + Gx1[92]*Gx2[90] + Gx1[106]*Gx2[104] + Gx1[120]*Gx2[118] + Gx1[134]*Gx2[132] + Gx1[148]*Gx2[146] + Gx1[162]*Gx2[160] + Gx1[176]*Gx2[174] + Gx1[190]*Gx2[188];
nmpcWorkspace.H[1719] += + Gx1[8]*Gx2[7] + Gx1[22]*Gx2[21] + Gx1[36]*Gx2[35] + Gx1[50]*Gx2[49] + Gx1[64]*Gx2[63] + Gx1[78]*Gx2[77] + Gx1[92]*Gx2[91] + Gx1[106]*Gx2[105] + Gx1[120]*Gx2[119] + Gx1[134]*Gx2[133] + Gx1[148]*Gx2[147] + Gx1[162]*Gx2[161] + Gx1[176]*Gx2[175] + Gx1[190]*Gx2[189];
nmpcWorkspace.H[1720] += + Gx1[8]*Gx2[8] + Gx1[22]*Gx2[22] + Gx1[36]*Gx2[36] + Gx1[50]*Gx2[50] + Gx1[64]*Gx2[64] + Gx1[78]*Gx2[78] + Gx1[92]*Gx2[92] + Gx1[106]*Gx2[106] + Gx1[120]*Gx2[120] + Gx1[134]*Gx2[134] + Gx1[148]*Gx2[148] + Gx1[162]*Gx2[162] + Gx1[176]*Gx2[176] + Gx1[190]*Gx2[190];
nmpcWorkspace.H[1721] += + Gx1[8]*Gx2[9] + Gx1[22]*Gx2[23] + Gx1[36]*Gx2[37] + Gx1[50]*Gx2[51] + Gx1[64]*Gx2[65] + Gx1[78]*Gx2[79] + Gx1[92]*Gx2[93] + Gx1[106]*Gx2[107] + Gx1[120]*Gx2[121] + Gx1[134]*Gx2[135] + Gx1[148]*Gx2[149] + Gx1[162]*Gx2[163] + Gx1[176]*Gx2[177] + Gx1[190]*Gx2[191];
nmpcWorkspace.H[1722] += + Gx1[8]*Gx2[10] + Gx1[22]*Gx2[24] + Gx1[36]*Gx2[38] + Gx1[50]*Gx2[52] + Gx1[64]*Gx2[66] + Gx1[78]*Gx2[80] + Gx1[92]*Gx2[94] + Gx1[106]*Gx2[108] + Gx1[120]*Gx2[122] + Gx1[134]*Gx2[136] + Gx1[148]*Gx2[150] + Gx1[162]*Gx2[164] + Gx1[176]*Gx2[178] + Gx1[190]*Gx2[192];
nmpcWorkspace.H[1723] += + Gx1[8]*Gx2[11] + Gx1[22]*Gx2[25] + Gx1[36]*Gx2[39] + Gx1[50]*Gx2[53] + Gx1[64]*Gx2[67] + Gx1[78]*Gx2[81] + Gx1[92]*Gx2[95] + Gx1[106]*Gx2[109] + Gx1[120]*Gx2[123] + Gx1[134]*Gx2[137] + Gx1[148]*Gx2[151] + Gx1[162]*Gx2[165] + Gx1[176]*Gx2[179] + Gx1[190]*Gx2[193];
nmpcWorkspace.H[1724] += + Gx1[8]*Gx2[12] + Gx1[22]*Gx2[26] + Gx1[36]*Gx2[40] + Gx1[50]*Gx2[54] + Gx1[64]*Gx2[68] + Gx1[78]*Gx2[82] + Gx1[92]*Gx2[96] + Gx1[106]*Gx2[110] + Gx1[120]*Gx2[124] + Gx1[134]*Gx2[138] + Gx1[148]*Gx2[152] + Gx1[162]*Gx2[166] + Gx1[176]*Gx2[180] + Gx1[190]*Gx2[194];
nmpcWorkspace.H[1725] += + Gx1[8]*Gx2[13] + Gx1[22]*Gx2[27] + Gx1[36]*Gx2[41] + Gx1[50]*Gx2[55] + Gx1[64]*Gx2[69] + Gx1[78]*Gx2[83] + Gx1[92]*Gx2[97] + Gx1[106]*Gx2[111] + Gx1[120]*Gx2[125] + Gx1[134]*Gx2[139] + Gx1[148]*Gx2[153] + Gx1[162]*Gx2[167] + Gx1[176]*Gx2[181] + Gx1[190]*Gx2[195];
nmpcWorkspace.H[1926] += + Gx1[9]*Gx2[0] + Gx1[23]*Gx2[14] + Gx1[37]*Gx2[28] + Gx1[51]*Gx2[42] + Gx1[65]*Gx2[56] + Gx1[79]*Gx2[70] + Gx1[93]*Gx2[84] + Gx1[107]*Gx2[98] + Gx1[121]*Gx2[112] + Gx1[135]*Gx2[126] + Gx1[149]*Gx2[140] + Gx1[163]*Gx2[154] + Gx1[177]*Gx2[168] + Gx1[191]*Gx2[182];
nmpcWorkspace.H[1927] += + Gx1[9]*Gx2[1] + Gx1[23]*Gx2[15] + Gx1[37]*Gx2[29] + Gx1[51]*Gx2[43] + Gx1[65]*Gx2[57] + Gx1[79]*Gx2[71] + Gx1[93]*Gx2[85] + Gx1[107]*Gx2[99] + Gx1[121]*Gx2[113] + Gx1[135]*Gx2[127] + Gx1[149]*Gx2[141] + Gx1[163]*Gx2[155] + Gx1[177]*Gx2[169] + Gx1[191]*Gx2[183];
nmpcWorkspace.H[1928] += + Gx1[9]*Gx2[2] + Gx1[23]*Gx2[16] + Gx1[37]*Gx2[30] + Gx1[51]*Gx2[44] + Gx1[65]*Gx2[58] + Gx1[79]*Gx2[72] + Gx1[93]*Gx2[86] + Gx1[107]*Gx2[100] + Gx1[121]*Gx2[114] + Gx1[135]*Gx2[128] + Gx1[149]*Gx2[142] + Gx1[163]*Gx2[156] + Gx1[177]*Gx2[170] + Gx1[191]*Gx2[184];
nmpcWorkspace.H[1929] += + Gx1[9]*Gx2[3] + Gx1[23]*Gx2[17] + Gx1[37]*Gx2[31] + Gx1[51]*Gx2[45] + Gx1[65]*Gx2[59] + Gx1[79]*Gx2[73] + Gx1[93]*Gx2[87] + Gx1[107]*Gx2[101] + Gx1[121]*Gx2[115] + Gx1[135]*Gx2[129] + Gx1[149]*Gx2[143] + Gx1[163]*Gx2[157] + Gx1[177]*Gx2[171] + Gx1[191]*Gx2[185];
nmpcWorkspace.H[1930] += + Gx1[9]*Gx2[4] + Gx1[23]*Gx2[18] + Gx1[37]*Gx2[32] + Gx1[51]*Gx2[46] + Gx1[65]*Gx2[60] + Gx1[79]*Gx2[74] + Gx1[93]*Gx2[88] + Gx1[107]*Gx2[102] + Gx1[121]*Gx2[116] + Gx1[135]*Gx2[130] + Gx1[149]*Gx2[144] + Gx1[163]*Gx2[158] + Gx1[177]*Gx2[172] + Gx1[191]*Gx2[186];
nmpcWorkspace.H[1931] += + Gx1[9]*Gx2[5] + Gx1[23]*Gx2[19] + Gx1[37]*Gx2[33] + Gx1[51]*Gx2[47] + Gx1[65]*Gx2[61] + Gx1[79]*Gx2[75] + Gx1[93]*Gx2[89] + Gx1[107]*Gx2[103] + Gx1[121]*Gx2[117] + Gx1[135]*Gx2[131] + Gx1[149]*Gx2[145] + Gx1[163]*Gx2[159] + Gx1[177]*Gx2[173] + Gx1[191]*Gx2[187];
nmpcWorkspace.H[1932] += + Gx1[9]*Gx2[6] + Gx1[23]*Gx2[20] + Gx1[37]*Gx2[34] + Gx1[51]*Gx2[48] + Gx1[65]*Gx2[62] + Gx1[79]*Gx2[76] + Gx1[93]*Gx2[90] + Gx1[107]*Gx2[104] + Gx1[121]*Gx2[118] + Gx1[135]*Gx2[132] + Gx1[149]*Gx2[146] + Gx1[163]*Gx2[160] + Gx1[177]*Gx2[174] + Gx1[191]*Gx2[188];
nmpcWorkspace.H[1933] += + Gx1[9]*Gx2[7] + Gx1[23]*Gx2[21] + Gx1[37]*Gx2[35] + Gx1[51]*Gx2[49] + Gx1[65]*Gx2[63] + Gx1[79]*Gx2[77] + Gx1[93]*Gx2[91] + Gx1[107]*Gx2[105] + Gx1[121]*Gx2[119] + Gx1[135]*Gx2[133] + Gx1[149]*Gx2[147] + Gx1[163]*Gx2[161] + Gx1[177]*Gx2[175] + Gx1[191]*Gx2[189];
nmpcWorkspace.H[1934] += + Gx1[9]*Gx2[8] + Gx1[23]*Gx2[22] + Gx1[37]*Gx2[36] + Gx1[51]*Gx2[50] + Gx1[65]*Gx2[64] + Gx1[79]*Gx2[78] + Gx1[93]*Gx2[92] + Gx1[107]*Gx2[106] + Gx1[121]*Gx2[120] + Gx1[135]*Gx2[134] + Gx1[149]*Gx2[148] + Gx1[163]*Gx2[162] + Gx1[177]*Gx2[176] + Gx1[191]*Gx2[190];
nmpcWorkspace.H[1935] += + Gx1[9]*Gx2[9] + Gx1[23]*Gx2[23] + Gx1[37]*Gx2[37] + Gx1[51]*Gx2[51] + Gx1[65]*Gx2[65] + Gx1[79]*Gx2[79] + Gx1[93]*Gx2[93] + Gx1[107]*Gx2[107] + Gx1[121]*Gx2[121] + Gx1[135]*Gx2[135] + Gx1[149]*Gx2[149] + Gx1[163]*Gx2[163] + Gx1[177]*Gx2[177] + Gx1[191]*Gx2[191];
nmpcWorkspace.H[1936] += + Gx1[9]*Gx2[10] + Gx1[23]*Gx2[24] + Gx1[37]*Gx2[38] + Gx1[51]*Gx2[52] + Gx1[65]*Gx2[66] + Gx1[79]*Gx2[80] + Gx1[93]*Gx2[94] + Gx1[107]*Gx2[108] + Gx1[121]*Gx2[122] + Gx1[135]*Gx2[136] + Gx1[149]*Gx2[150] + Gx1[163]*Gx2[164] + Gx1[177]*Gx2[178] + Gx1[191]*Gx2[192];
nmpcWorkspace.H[1937] += + Gx1[9]*Gx2[11] + Gx1[23]*Gx2[25] + Gx1[37]*Gx2[39] + Gx1[51]*Gx2[53] + Gx1[65]*Gx2[67] + Gx1[79]*Gx2[81] + Gx1[93]*Gx2[95] + Gx1[107]*Gx2[109] + Gx1[121]*Gx2[123] + Gx1[135]*Gx2[137] + Gx1[149]*Gx2[151] + Gx1[163]*Gx2[165] + Gx1[177]*Gx2[179] + Gx1[191]*Gx2[193];
nmpcWorkspace.H[1938] += + Gx1[9]*Gx2[12] + Gx1[23]*Gx2[26] + Gx1[37]*Gx2[40] + Gx1[51]*Gx2[54] + Gx1[65]*Gx2[68] + Gx1[79]*Gx2[82] + Gx1[93]*Gx2[96] + Gx1[107]*Gx2[110] + Gx1[121]*Gx2[124] + Gx1[135]*Gx2[138] + Gx1[149]*Gx2[152] + Gx1[163]*Gx2[166] + Gx1[177]*Gx2[180] + Gx1[191]*Gx2[194];
nmpcWorkspace.H[1939] += + Gx1[9]*Gx2[13] + Gx1[23]*Gx2[27] + Gx1[37]*Gx2[41] + Gx1[51]*Gx2[55] + Gx1[65]*Gx2[69] + Gx1[79]*Gx2[83] + Gx1[93]*Gx2[97] + Gx1[107]*Gx2[111] + Gx1[121]*Gx2[125] + Gx1[135]*Gx2[139] + Gx1[149]*Gx2[153] + Gx1[163]*Gx2[167] + Gx1[177]*Gx2[181] + Gx1[191]*Gx2[195];
nmpcWorkspace.H[2140] += + Gx1[10]*Gx2[0] + Gx1[24]*Gx2[14] + Gx1[38]*Gx2[28] + Gx1[52]*Gx2[42] + Gx1[66]*Gx2[56] + Gx1[80]*Gx2[70] + Gx1[94]*Gx2[84] + Gx1[108]*Gx2[98] + Gx1[122]*Gx2[112] + Gx1[136]*Gx2[126] + Gx1[150]*Gx2[140] + Gx1[164]*Gx2[154] + Gx1[178]*Gx2[168] + Gx1[192]*Gx2[182];
nmpcWorkspace.H[2141] += + Gx1[10]*Gx2[1] + Gx1[24]*Gx2[15] + Gx1[38]*Gx2[29] + Gx1[52]*Gx2[43] + Gx1[66]*Gx2[57] + Gx1[80]*Gx2[71] + Gx1[94]*Gx2[85] + Gx1[108]*Gx2[99] + Gx1[122]*Gx2[113] + Gx1[136]*Gx2[127] + Gx1[150]*Gx2[141] + Gx1[164]*Gx2[155] + Gx1[178]*Gx2[169] + Gx1[192]*Gx2[183];
nmpcWorkspace.H[2142] += + Gx1[10]*Gx2[2] + Gx1[24]*Gx2[16] + Gx1[38]*Gx2[30] + Gx1[52]*Gx2[44] + Gx1[66]*Gx2[58] + Gx1[80]*Gx2[72] + Gx1[94]*Gx2[86] + Gx1[108]*Gx2[100] + Gx1[122]*Gx2[114] + Gx1[136]*Gx2[128] + Gx1[150]*Gx2[142] + Gx1[164]*Gx2[156] + Gx1[178]*Gx2[170] + Gx1[192]*Gx2[184];
nmpcWorkspace.H[2143] += + Gx1[10]*Gx2[3] + Gx1[24]*Gx2[17] + Gx1[38]*Gx2[31] + Gx1[52]*Gx2[45] + Gx1[66]*Gx2[59] + Gx1[80]*Gx2[73] + Gx1[94]*Gx2[87] + Gx1[108]*Gx2[101] + Gx1[122]*Gx2[115] + Gx1[136]*Gx2[129] + Gx1[150]*Gx2[143] + Gx1[164]*Gx2[157] + Gx1[178]*Gx2[171] + Gx1[192]*Gx2[185];
nmpcWorkspace.H[2144] += + Gx1[10]*Gx2[4] + Gx1[24]*Gx2[18] + Gx1[38]*Gx2[32] + Gx1[52]*Gx2[46] + Gx1[66]*Gx2[60] + Gx1[80]*Gx2[74] + Gx1[94]*Gx2[88] + Gx1[108]*Gx2[102] + Gx1[122]*Gx2[116] + Gx1[136]*Gx2[130] + Gx1[150]*Gx2[144] + Gx1[164]*Gx2[158] + Gx1[178]*Gx2[172] + Gx1[192]*Gx2[186];
nmpcWorkspace.H[2145] += + Gx1[10]*Gx2[5] + Gx1[24]*Gx2[19] + Gx1[38]*Gx2[33] + Gx1[52]*Gx2[47] + Gx1[66]*Gx2[61] + Gx1[80]*Gx2[75] + Gx1[94]*Gx2[89] + Gx1[108]*Gx2[103] + Gx1[122]*Gx2[117] + Gx1[136]*Gx2[131] + Gx1[150]*Gx2[145] + Gx1[164]*Gx2[159] + Gx1[178]*Gx2[173] + Gx1[192]*Gx2[187];
nmpcWorkspace.H[2146] += + Gx1[10]*Gx2[6] + Gx1[24]*Gx2[20] + Gx1[38]*Gx2[34] + Gx1[52]*Gx2[48] + Gx1[66]*Gx2[62] + Gx1[80]*Gx2[76] + Gx1[94]*Gx2[90] + Gx1[108]*Gx2[104] + Gx1[122]*Gx2[118] + Gx1[136]*Gx2[132] + Gx1[150]*Gx2[146] + Gx1[164]*Gx2[160] + Gx1[178]*Gx2[174] + Gx1[192]*Gx2[188];
nmpcWorkspace.H[2147] += + Gx1[10]*Gx2[7] + Gx1[24]*Gx2[21] + Gx1[38]*Gx2[35] + Gx1[52]*Gx2[49] + Gx1[66]*Gx2[63] + Gx1[80]*Gx2[77] + Gx1[94]*Gx2[91] + Gx1[108]*Gx2[105] + Gx1[122]*Gx2[119] + Gx1[136]*Gx2[133] + Gx1[150]*Gx2[147] + Gx1[164]*Gx2[161] + Gx1[178]*Gx2[175] + Gx1[192]*Gx2[189];
nmpcWorkspace.H[2148] += + Gx1[10]*Gx2[8] + Gx1[24]*Gx2[22] + Gx1[38]*Gx2[36] + Gx1[52]*Gx2[50] + Gx1[66]*Gx2[64] + Gx1[80]*Gx2[78] + Gx1[94]*Gx2[92] + Gx1[108]*Gx2[106] + Gx1[122]*Gx2[120] + Gx1[136]*Gx2[134] + Gx1[150]*Gx2[148] + Gx1[164]*Gx2[162] + Gx1[178]*Gx2[176] + Gx1[192]*Gx2[190];
nmpcWorkspace.H[2149] += + Gx1[10]*Gx2[9] + Gx1[24]*Gx2[23] + Gx1[38]*Gx2[37] + Gx1[52]*Gx2[51] + Gx1[66]*Gx2[65] + Gx1[80]*Gx2[79] + Gx1[94]*Gx2[93] + Gx1[108]*Gx2[107] + Gx1[122]*Gx2[121] + Gx1[136]*Gx2[135] + Gx1[150]*Gx2[149] + Gx1[164]*Gx2[163] + Gx1[178]*Gx2[177] + Gx1[192]*Gx2[191];
nmpcWorkspace.H[2150] += + Gx1[10]*Gx2[10] + Gx1[24]*Gx2[24] + Gx1[38]*Gx2[38] + Gx1[52]*Gx2[52] + Gx1[66]*Gx2[66] + Gx1[80]*Gx2[80] + Gx1[94]*Gx2[94] + Gx1[108]*Gx2[108] + Gx1[122]*Gx2[122] + Gx1[136]*Gx2[136] + Gx1[150]*Gx2[150] + Gx1[164]*Gx2[164] + Gx1[178]*Gx2[178] + Gx1[192]*Gx2[192];
nmpcWorkspace.H[2151] += + Gx1[10]*Gx2[11] + Gx1[24]*Gx2[25] + Gx1[38]*Gx2[39] + Gx1[52]*Gx2[53] + Gx1[66]*Gx2[67] + Gx1[80]*Gx2[81] + Gx1[94]*Gx2[95] + Gx1[108]*Gx2[109] + Gx1[122]*Gx2[123] + Gx1[136]*Gx2[137] + Gx1[150]*Gx2[151] + Gx1[164]*Gx2[165] + Gx1[178]*Gx2[179] + Gx1[192]*Gx2[193];
nmpcWorkspace.H[2152] += + Gx1[10]*Gx2[12] + Gx1[24]*Gx2[26] + Gx1[38]*Gx2[40] + Gx1[52]*Gx2[54] + Gx1[66]*Gx2[68] + Gx1[80]*Gx2[82] + Gx1[94]*Gx2[96] + Gx1[108]*Gx2[110] + Gx1[122]*Gx2[124] + Gx1[136]*Gx2[138] + Gx1[150]*Gx2[152] + Gx1[164]*Gx2[166] + Gx1[178]*Gx2[180] + Gx1[192]*Gx2[194];
nmpcWorkspace.H[2153] += + Gx1[10]*Gx2[13] + Gx1[24]*Gx2[27] + Gx1[38]*Gx2[41] + Gx1[52]*Gx2[55] + Gx1[66]*Gx2[69] + Gx1[80]*Gx2[83] + Gx1[94]*Gx2[97] + Gx1[108]*Gx2[111] + Gx1[122]*Gx2[125] + Gx1[136]*Gx2[139] + Gx1[150]*Gx2[153] + Gx1[164]*Gx2[167] + Gx1[178]*Gx2[181] + Gx1[192]*Gx2[195];
nmpcWorkspace.H[2354] += + Gx1[11]*Gx2[0] + Gx1[25]*Gx2[14] + Gx1[39]*Gx2[28] + Gx1[53]*Gx2[42] + Gx1[67]*Gx2[56] + Gx1[81]*Gx2[70] + Gx1[95]*Gx2[84] + Gx1[109]*Gx2[98] + Gx1[123]*Gx2[112] + Gx1[137]*Gx2[126] + Gx1[151]*Gx2[140] + Gx1[165]*Gx2[154] + Gx1[179]*Gx2[168] + Gx1[193]*Gx2[182];
nmpcWorkspace.H[2355] += + Gx1[11]*Gx2[1] + Gx1[25]*Gx2[15] + Gx1[39]*Gx2[29] + Gx1[53]*Gx2[43] + Gx1[67]*Gx2[57] + Gx1[81]*Gx2[71] + Gx1[95]*Gx2[85] + Gx1[109]*Gx2[99] + Gx1[123]*Gx2[113] + Gx1[137]*Gx2[127] + Gx1[151]*Gx2[141] + Gx1[165]*Gx2[155] + Gx1[179]*Gx2[169] + Gx1[193]*Gx2[183];
nmpcWorkspace.H[2356] += + Gx1[11]*Gx2[2] + Gx1[25]*Gx2[16] + Gx1[39]*Gx2[30] + Gx1[53]*Gx2[44] + Gx1[67]*Gx2[58] + Gx1[81]*Gx2[72] + Gx1[95]*Gx2[86] + Gx1[109]*Gx2[100] + Gx1[123]*Gx2[114] + Gx1[137]*Gx2[128] + Gx1[151]*Gx2[142] + Gx1[165]*Gx2[156] + Gx1[179]*Gx2[170] + Gx1[193]*Gx2[184];
nmpcWorkspace.H[2357] += + Gx1[11]*Gx2[3] + Gx1[25]*Gx2[17] + Gx1[39]*Gx2[31] + Gx1[53]*Gx2[45] + Gx1[67]*Gx2[59] + Gx1[81]*Gx2[73] + Gx1[95]*Gx2[87] + Gx1[109]*Gx2[101] + Gx1[123]*Gx2[115] + Gx1[137]*Gx2[129] + Gx1[151]*Gx2[143] + Gx1[165]*Gx2[157] + Gx1[179]*Gx2[171] + Gx1[193]*Gx2[185];
nmpcWorkspace.H[2358] += + Gx1[11]*Gx2[4] + Gx1[25]*Gx2[18] + Gx1[39]*Gx2[32] + Gx1[53]*Gx2[46] + Gx1[67]*Gx2[60] + Gx1[81]*Gx2[74] + Gx1[95]*Gx2[88] + Gx1[109]*Gx2[102] + Gx1[123]*Gx2[116] + Gx1[137]*Gx2[130] + Gx1[151]*Gx2[144] + Gx1[165]*Gx2[158] + Gx1[179]*Gx2[172] + Gx1[193]*Gx2[186];
nmpcWorkspace.H[2359] += + Gx1[11]*Gx2[5] + Gx1[25]*Gx2[19] + Gx1[39]*Gx2[33] + Gx1[53]*Gx2[47] + Gx1[67]*Gx2[61] + Gx1[81]*Gx2[75] + Gx1[95]*Gx2[89] + Gx1[109]*Gx2[103] + Gx1[123]*Gx2[117] + Gx1[137]*Gx2[131] + Gx1[151]*Gx2[145] + Gx1[165]*Gx2[159] + Gx1[179]*Gx2[173] + Gx1[193]*Gx2[187];
nmpcWorkspace.H[2360] += + Gx1[11]*Gx2[6] + Gx1[25]*Gx2[20] + Gx1[39]*Gx2[34] + Gx1[53]*Gx2[48] + Gx1[67]*Gx2[62] + Gx1[81]*Gx2[76] + Gx1[95]*Gx2[90] + Gx1[109]*Gx2[104] + Gx1[123]*Gx2[118] + Gx1[137]*Gx2[132] + Gx1[151]*Gx2[146] + Gx1[165]*Gx2[160] + Gx1[179]*Gx2[174] + Gx1[193]*Gx2[188];
nmpcWorkspace.H[2361] += + Gx1[11]*Gx2[7] + Gx1[25]*Gx2[21] + Gx1[39]*Gx2[35] + Gx1[53]*Gx2[49] + Gx1[67]*Gx2[63] + Gx1[81]*Gx2[77] + Gx1[95]*Gx2[91] + Gx1[109]*Gx2[105] + Gx1[123]*Gx2[119] + Gx1[137]*Gx2[133] + Gx1[151]*Gx2[147] + Gx1[165]*Gx2[161] + Gx1[179]*Gx2[175] + Gx1[193]*Gx2[189];
nmpcWorkspace.H[2362] += + Gx1[11]*Gx2[8] + Gx1[25]*Gx2[22] + Gx1[39]*Gx2[36] + Gx1[53]*Gx2[50] + Gx1[67]*Gx2[64] + Gx1[81]*Gx2[78] + Gx1[95]*Gx2[92] + Gx1[109]*Gx2[106] + Gx1[123]*Gx2[120] + Gx1[137]*Gx2[134] + Gx1[151]*Gx2[148] + Gx1[165]*Gx2[162] + Gx1[179]*Gx2[176] + Gx1[193]*Gx2[190];
nmpcWorkspace.H[2363] += + Gx1[11]*Gx2[9] + Gx1[25]*Gx2[23] + Gx1[39]*Gx2[37] + Gx1[53]*Gx2[51] + Gx1[67]*Gx2[65] + Gx1[81]*Gx2[79] + Gx1[95]*Gx2[93] + Gx1[109]*Gx2[107] + Gx1[123]*Gx2[121] + Gx1[137]*Gx2[135] + Gx1[151]*Gx2[149] + Gx1[165]*Gx2[163] + Gx1[179]*Gx2[177] + Gx1[193]*Gx2[191];
nmpcWorkspace.H[2364] += + Gx1[11]*Gx2[10] + Gx1[25]*Gx2[24] + Gx1[39]*Gx2[38] + Gx1[53]*Gx2[52] + Gx1[67]*Gx2[66] + Gx1[81]*Gx2[80] + Gx1[95]*Gx2[94] + Gx1[109]*Gx2[108] + Gx1[123]*Gx2[122] + Gx1[137]*Gx2[136] + Gx1[151]*Gx2[150] + Gx1[165]*Gx2[164] + Gx1[179]*Gx2[178] + Gx1[193]*Gx2[192];
nmpcWorkspace.H[2365] += + Gx1[11]*Gx2[11] + Gx1[25]*Gx2[25] + Gx1[39]*Gx2[39] + Gx1[53]*Gx2[53] + Gx1[67]*Gx2[67] + Gx1[81]*Gx2[81] + Gx1[95]*Gx2[95] + Gx1[109]*Gx2[109] + Gx1[123]*Gx2[123] + Gx1[137]*Gx2[137] + Gx1[151]*Gx2[151] + Gx1[165]*Gx2[165] + Gx1[179]*Gx2[179] + Gx1[193]*Gx2[193];
nmpcWorkspace.H[2366] += + Gx1[11]*Gx2[12] + Gx1[25]*Gx2[26] + Gx1[39]*Gx2[40] + Gx1[53]*Gx2[54] + Gx1[67]*Gx2[68] + Gx1[81]*Gx2[82] + Gx1[95]*Gx2[96] + Gx1[109]*Gx2[110] + Gx1[123]*Gx2[124] + Gx1[137]*Gx2[138] + Gx1[151]*Gx2[152] + Gx1[165]*Gx2[166] + Gx1[179]*Gx2[180] + Gx1[193]*Gx2[194];
nmpcWorkspace.H[2367] += + Gx1[11]*Gx2[13] + Gx1[25]*Gx2[27] + Gx1[39]*Gx2[41] + Gx1[53]*Gx2[55] + Gx1[67]*Gx2[69] + Gx1[81]*Gx2[83] + Gx1[95]*Gx2[97] + Gx1[109]*Gx2[111] + Gx1[123]*Gx2[125] + Gx1[137]*Gx2[139] + Gx1[151]*Gx2[153] + Gx1[165]*Gx2[167] + Gx1[179]*Gx2[181] + Gx1[193]*Gx2[195];
nmpcWorkspace.H[2568] += + Gx1[12]*Gx2[0] + Gx1[26]*Gx2[14] + Gx1[40]*Gx2[28] + Gx1[54]*Gx2[42] + Gx1[68]*Gx2[56] + Gx1[82]*Gx2[70] + Gx1[96]*Gx2[84] + Gx1[110]*Gx2[98] + Gx1[124]*Gx2[112] + Gx1[138]*Gx2[126] + Gx1[152]*Gx2[140] + Gx1[166]*Gx2[154] + Gx1[180]*Gx2[168] + Gx1[194]*Gx2[182];
nmpcWorkspace.H[2569] += + Gx1[12]*Gx2[1] + Gx1[26]*Gx2[15] + Gx1[40]*Gx2[29] + Gx1[54]*Gx2[43] + Gx1[68]*Gx2[57] + Gx1[82]*Gx2[71] + Gx1[96]*Gx2[85] + Gx1[110]*Gx2[99] + Gx1[124]*Gx2[113] + Gx1[138]*Gx2[127] + Gx1[152]*Gx2[141] + Gx1[166]*Gx2[155] + Gx1[180]*Gx2[169] + Gx1[194]*Gx2[183];
nmpcWorkspace.H[2570] += + Gx1[12]*Gx2[2] + Gx1[26]*Gx2[16] + Gx1[40]*Gx2[30] + Gx1[54]*Gx2[44] + Gx1[68]*Gx2[58] + Gx1[82]*Gx2[72] + Gx1[96]*Gx2[86] + Gx1[110]*Gx2[100] + Gx1[124]*Gx2[114] + Gx1[138]*Gx2[128] + Gx1[152]*Gx2[142] + Gx1[166]*Gx2[156] + Gx1[180]*Gx2[170] + Gx1[194]*Gx2[184];
nmpcWorkspace.H[2571] += + Gx1[12]*Gx2[3] + Gx1[26]*Gx2[17] + Gx1[40]*Gx2[31] + Gx1[54]*Gx2[45] + Gx1[68]*Gx2[59] + Gx1[82]*Gx2[73] + Gx1[96]*Gx2[87] + Gx1[110]*Gx2[101] + Gx1[124]*Gx2[115] + Gx1[138]*Gx2[129] + Gx1[152]*Gx2[143] + Gx1[166]*Gx2[157] + Gx1[180]*Gx2[171] + Gx1[194]*Gx2[185];
nmpcWorkspace.H[2572] += + Gx1[12]*Gx2[4] + Gx1[26]*Gx2[18] + Gx1[40]*Gx2[32] + Gx1[54]*Gx2[46] + Gx1[68]*Gx2[60] + Gx1[82]*Gx2[74] + Gx1[96]*Gx2[88] + Gx1[110]*Gx2[102] + Gx1[124]*Gx2[116] + Gx1[138]*Gx2[130] + Gx1[152]*Gx2[144] + Gx1[166]*Gx2[158] + Gx1[180]*Gx2[172] + Gx1[194]*Gx2[186];
nmpcWorkspace.H[2573] += + Gx1[12]*Gx2[5] + Gx1[26]*Gx2[19] + Gx1[40]*Gx2[33] + Gx1[54]*Gx2[47] + Gx1[68]*Gx2[61] + Gx1[82]*Gx2[75] + Gx1[96]*Gx2[89] + Gx1[110]*Gx2[103] + Gx1[124]*Gx2[117] + Gx1[138]*Gx2[131] + Gx1[152]*Gx2[145] + Gx1[166]*Gx2[159] + Gx1[180]*Gx2[173] + Gx1[194]*Gx2[187];
nmpcWorkspace.H[2574] += + Gx1[12]*Gx2[6] + Gx1[26]*Gx2[20] + Gx1[40]*Gx2[34] + Gx1[54]*Gx2[48] + Gx1[68]*Gx2[62] + Gx1[82]*Gx2[76] + Gx1[96]*Gx2[90] + Gx1[110]*Gx2[104] + Gx1[124]*Gx2[118] + Gx1[138]*Gx2[132] + Gx1[152]*Gx2[146] + Gx1[166]*Gx2[160] + Gx1[180]*Gx2[174] + Gx1[194]*Gx2[188];
nmpcWorkspace.H[2575] += + Gx1[12]*Gx2[7] + Gx1[26]*Gx2[21] + Gx1[40]*Gx2[35] + Gx1[54]*Gx2[49] + Gx1[68]*Gx2[63] + Gx1[82]*Gx2[77] + Gx1[96]*Gx2[91] + Gx1[110]*Gx2[105] + Gx1[124]*Gx2[119] + Gx1[138]*Gx2[133] + Gx1[152]*Gx2[147] + Gx1[166]*Gx2[161] + Gx1[180]*Gx2[175] + Gx1[194]*Gx2[189];
nmpcWorkspace.H[2576] += + Gx1[12]*Gx2[8] + Gx1[26]*Gx2[22] + Gx1[40]*Gx2[36] + Gx1[54]*Gx2[50] + Gx1[68]*Gx2[64] + Gx1[82]*Gx2[78] + Gx1[96]*Gx2[92] + Gx1[110]*Gx2[106] + Gx1[124]*Gx2[120] + Gx1[138]*Gx2[134] + Gx1[152]*Gx2[148] + Gx1[166]*Gx2[162] + Gx1[180]*Gx2[176] + Gx1[194]*Gx2[190];
nmpcWorkspace.H[2577] += + Gx1[12]*Gx2[9] + Gx1[26]*Gx2[23] + Gx1[40]*Gx2[37] + Gx1[54]*Gx2[51] + Gx1[68]*Gx2[65] + Gx1[82]*Gx2[79] + Gx1[96]*Gx2[93] + Gx1[110]*Gx2[107] + Gx1[124]*Gx2[121] + Gx1[138]*Gx2[135] + Gx1[152]*Gx2[149] + Gx1[166]*Gx2[163] + Gx1[180]*Gx2[177] + Gx1[194]*Gx2[191];
nmpcWorkspace.H[2578] += + Gx1[12]*Gx2[10] + Gx1[26]*Gx2[24] + Gx1[40]*Gx2[38] + Gx1[54]*Gx2[52] + Gx1[68]*Gx2[66] + Gx1[82]*Gx2[80] + Gx1[96]*Gx2[94] + Gx1[110]*Gx2[108] + Gx1[124]*Gx2[122] + Gx1[138]*Gx2[136] + Gx1[152]*Gx2[150] + Gx1[166]*Gx2[164] + Gx1[180]*Gx2[178] + Gx1[194]*Gx2[192];
nmpcWorkspace.H[2579] += + Gx1[12]*Gx2[11] + Gx1[26]*Gx2[25] + Gx1[40]*Gx2[39] + Gx1[54]*Gx2[53] + Gx1[68]*Gx2[67] + Gx1[82]*Gx2[81] + Gx1[96]*Gx2[95] + Gx1[110]*Gx2[109] + Gx1[124]*Gx2[123] + Gx1[138]*Gx2[137] + Gx1[152]*Gx2[151] + Gx1[166]*Gx2[165] + Gx1[180]*Gx2[179] + Gx1[194]*Gx2[193];
nmpcWorkspace.H[2580] += + Gx1[12]*Gx2[12] + Gx1[26]*Gx2[26] + Gx1[40]*Gx2[40] + Gx1[54]*Gx2[54] + Gx1[68]*Gx2[68] + Gx1[82]*Gx2[82] + Gx1[96]*Gx2[96] + Gx1[110]*Gx2[110] + Gx1[124]*Gx2[124] + Gx1[138]*Gx2[138] + Gx1[152]*Gx2[152] + Gx1[166]*Gx2[166] + Gx1[180]*Gx2[180] + Gx1[194]*Gx2[194];
nmpcWorkspace.H[2581] += + Gx1[12]*Gx2[13] + Gx1[26]*Gx2[27] + Gx1[40]*Gx2[41] + Gx1[54]*Gx2[55] + Gx1[68]*Gx2[69] + Gx1[82]*Gx2[83] + Gx1[96]*Gx2[97] + Gx1[110]*Gx2[111] + Gx1[124]*Gx2[125] + Gx1[138]*Gx2[139] + Gx1[152]*Gx2[153] + Gx1[166]*Gx2[167] + Gx1[180]*Gx2[181] + Gx1[194]*Gx2[195];
nmpcWorkspace.H[2782] += + Gx1[13]*Gx2[0] + Gx1[27]*Gx2[14] + Gx1[41]*Gx2[28] + Gx1[55]*Gx2[42] + Gx1[69]*Gx2[56] + Gx1[83]*Gx2[70] + Gx1[97]*Gx2[84] + Gx1[111]*Gx2[98] + Gx1[125]*Gx2[112] + Gx1[139]*Gx2[126] + Gx1[153]*Gx2[140] + Gx1[167]*Gx2[154] + Gx1[181]*Gx2[168] + Gx1[195]*Gx2[182];
nmpcWorkspace.H[2783] += + Gx1[13]*Gx2[1] + Gx1[27]*Gx2[15] + Gx1[41]*Gx2[29] + Gx1[55]*Gx2[43] + Gx1[69]*Gx2[57] + Gx1[83]*Gx2[71] + Gx1[97]*Gx2[85] + Gx1[111]*Gx2[99] + Gx1[125]*Gx2[113] + Gx1[139]*Gx2[127] + Gx1[153]*Gx2[141] + Gx1[167]*Gx2[155] + Gx1[181]*Gx2[169] + Gx1[195]*Gx2[183];
nmpcWorkspace.H[2784] += + Gx1[13]*Gx2[2] + Gx1[27]*Gx2[16] + Gx1[41]*Gx2[30] + Gx1[55]*Gx2[44] + Gx1[69]*Gx2[58] + Gx1[83]*Gx2[72] + Gx1[97]*Gx2[86] + Gx1[111]*Gx2[100] + Gx1[125]*Gx2[114] + Gx1[139]*Gx2[128] + Gx1[153]*Gx2[142] + Gx1[167]*Gx2[156] + Gx1[181]*Gx2[170] + Gx1[195]*Gx2[184];
nmpcWorkspace.H[2785] += + Gx1[13]*Gx2[3] + Gx1[27]*Gx2[17] + Gx1[41]*Gx2[31] + Gx1[55]*Gx2[45] + Gx1[69]*Gx2[59] + Gx1[83]*Gx2[73] + Gx1[97]*Gx2[87] + Gx1[111]*Gx2[101] + Gx1[125]*Gx2[115] + Gx1[139]*Gx2[129] + Gx1[153]*Gx2[143] + Gx1[167]*Gx2[157] + Gx1[181]*Gx2[171] + Gx1[195]*Gx2[185];
nmpcWorkspace.H[2786] += + Gx1[13]*Gx2[4] + Gx1[27]*Gx2[18] + Gx1[41]*Gx2[32] + Gx1[55]*Gx2[46] + Gx1[69]*Gx2[60] + Gx1[83]*Gx2[74] + Gx1[97]*Gx2[88] + Gx1[111]*Gx2[102] + Gx1[125]*Gx2[116] + Gx1[139]*Gx2[130] + Gx1[153]*Gx2[144] + Gx1[167]*Gx2[158] + Gx1[181]*Gx2[172] + Gx1[195]*Gx2[186];
nmpcWorkspace.H[2787] += + Gx1[13]*Gx2[5] + Gx1[27]*Gx2[19] + Gx1[41]*Gx2[33] + Gx1[55]*Gx2[47] + Gx1[69]*Gx2[61] + Gx1[83]*Gx2[75] + Gx1[97]*Gx2[89] + Gx1[111]*Gx2[103] + Gx1[125]*Gx2[117] + Gx1[139]*Gx2[131] + Gx1[153]*Gx2[145] + Gx1[167]*Gx2[159] + Gx1[181]*Gx2[173] + Gx1[195]*Gx2[187];
nmpcWorkspace.H[2788] += + Gx1[13]*Gx2[6] + Gx1[27]*Gx2[20] + Gx1[41]*Gx2[34] + Gx1[55]*Gx2[48] + Gx1[69]*Gx2[62] + Gx1[83]*Gx2[76] + Gx1[97]*Gx2[90] + Gx1[111]*Gx2[104] + Gx1[125]*Gx2[118] + Gx1[139]*Gx2[132] + Gx1[153]*Gx2[146] + Gx1[167]*Gx2[160] + Gx1[181]*Gx2[174] + Gx1[195]*Gx2[188];
nmpcWorkspace.H[2789] += + Gx1[13]*Gx2[7] + Gx1[27]*Gx2[21] + Gx1[41]*Gx2[35] + Gx1[55]*Gx2[49] + Gx1[69]*Gx2[63] + Gx1[83]*Gx2[77] + Gx1[97]*Gx2[91] + Gx1[111]*Gx2[105] + Gx1[125]*Gx2[119] + Gx1[139]*Gx2[133] + Gx1[153]*Gx2[147] + Gx1[167]*Gx2[161] + Gx1[181]*Gx2[175] + Gx1[195]*Gx2[189];
nmpcWorkspace.H[2790] += + Gx1[13]*Gx2[8] + Gx1[27]*Gx2[22] + Gx1[41]*Gx2[36] + Gx1[55]*Gx2[50] + Gx1[69]*Gx2[64] + Gx1[83]*Gx2[78] + Gx1[97]*Gx2[92] + Gx1[111]*Gx2[106] + Gx1[125]*Gx2[120] + Gx1[139]*Gx2[134] + Gx1[153]*Gx2[148] + Gx1[167]*Gx2[162] + Gx1[181]*Gx2[176] + Gx1[195]*Gx2[190];
nmpcWorkspace.H[2791] += + Gx1[13]*Gx2[9] + Gx1[27]*Gx2[23] + Gx1[41]*Gx2[37] + Gx1[55]*Gx2[51] + Gx1[69]*Gx2[65] + Gx1[83]*Gx2[79] + Gx1[97]*Gx2[93] + Gx1[111]*Gx2[107] + Gx1[125]*Gx2[121] + Gx1[139]*Gx2[135] + Gx1[153]*Gx2[149] + Gx1[167]*Gx2[163] + Gx1[181]*Gx2[177] + Gx1[195]*Gx2[191];
nmpcWorkspace.H[2792] += + Gx1[13]*Gx2[10] + Gx1[27]*Gx2[24] + Gx1[41]*Gx2[38] + Gx1[55]*Gx2[52] + Gx1[69]*Gx2[66] + Gx1[83]*Gx2[80] + Gx1[97]*Gx2[94] + Gx1[111]*Gx2[108] + Gx1[125]*Gx2[122] + Gx1[139]*Gx2[136] + Gx1[153]*Gx2[150] + Gx1[167]*Gx2[164] + Gx1[181]*Gx2[178] + Gx1[195]*Gx2[192];
nmpcWorkspace.H[2793] += + Gx1[13]*Gx2[11] + Gx1[27]*Gx2[25] + Gx1[41]*Gx2[39] + Gx1[55]*Gx2[53] + Gx1[69]*Gx2[67] + Gx1[83]*Gx2[81] + Gx1[97]*Gx2[95] + Gx1[111]*Gx2[109] + Gx1[125]*Gx2[123] + Gx1[139]*Gx2[137] + Gx1[153]*Gx2[151] + Gx1[167]*Gx2[165] + Gx1[181]*Gx2[179] + Gx1[195]*Gx2[193];
nmpcWorkspace.H[2794] += + Gx1[13]*Gx2[12] + Gx1[27]*Gx2[26] + Gx1[41]*Gx2[40] + Gx1[55]*Gx2[54] + Gx1[69]*Gx2[68] + Gx1[83]*Gx2[82] + Gx1[97]*Gx2[96] + Gx1[111]*Gx2[110] + Gx1[125]*Gx2[124] + Gx1[139]*Gx2[138] + Gx1[153]*Gx2[152] + Gx1[167]*Gx2[166] + Gx1[181]*Gx2[180] + Gx1[195]*Gx2[194];
nmpcWorkspace.H[2795] += + Gx1[13]*Gx2[13] + Gx1[27]*Gx2[27] + Gx1[41]*Gx2[41] + Gx1[55]*Gx2[55] + Gx1[69]*Gx2[69] + Gx1[83]*Gx2[83] + Gx1[97]*Gx2[97] + Gx1[111]*Gx2[111] + Gx1[125]*Gx2[125] + Gx1[139]*Gx2[139] + Gx1[153]*Gx2[153] + Gx1[167]*Gx2[167] + Gx1[181]*Gx2[181] + Gx1[195]*Gx2[195];
}

void nmpc_macCTSlx( real_t* const C0, real_t* const g0 )
{
g0[0] += 0.0;
;
g0[1] += 0.0;
;
g0[2] += 0.0;
;
g0[3] += 0.0;
;
g0[4] += 0.0;
;
g0[5] += 0.0;
;
g0[6] += 0.0;
;
g0[7] += 0.0;
;
g0[8] += 0.0;
;
g0[9] += 0.0;
;
g0[10] += 0.0;
;
g0[11] += 0.0;
;
g0[12] += 0.0;
;
g0[13] += 0.0;
;
}

void nmpc_macETSlu( real_t* const E0, real_t* const g1 )
{
g1[0] += 0.0;
;
g1[1] += 0.0;
;
g1[2] += 0.0;
;
g1[3] += 0.0;
;
}

void nmpc_condensePrep(  )
{
int lRun1;
int lRun2;
int lRun3;
int lRun4;
int lRun5;
nmpc_moveGuE( nmpcWorkspace.evGu, nmpcWorkspace.E );
for (lRun1 = 1; lRun1 < 50; ++lRun1)
{
nmpc_moveGxT( &(nmpcWorkspace.evGx[ lRun1 * 196 ]), nmpcWorkspace.T );
nmpc_multGxd( &(nmpcWorkspace.d[ lRun1 * 14-14 ]), &(nmpcWorkspace.evGx[ lRun1 * 196 ]), &(nmpcWorkspace.d[ lRun1 * 14 ]) );
nmpc_multGxGx( nmpcWorkspace.T, &(nmpcWorkspace.evGx[ lRun1 * 196-196 ]), &(nmpcWorkspace.evGx[ lRun1 * 196 ]) );
for (lRun2 = 0; lRun2 < lRun1; ++lRun2)
{
lRun4 = (((lRun1) * (lRun1-1)) / (2)) + (lRun2);
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
nmpc_multGxGu( nmpcWorkspace.T, &(nmpcWorkspace.E[ lRun4 * 56 ]), &(nmpcWorkspace.E[ lRun3 * 56 ]) );
}
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
nmpc_moveGuE( &(nmpcWorkspace.evGu[ lRun1 * 56 ]), &(nmpcWorkspace.E[ lRun3 * 56 ]) );
}

nmpc_multGxGx( &(nmpcWorkspace.Q1[ 196 ]), nmpcWorkspace.evGx, nmpcWorkspace.QGx );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 392 ]), &(nmpcWorkspace.evGx[ 196 ]), &(nmpcWorkspace.QGx[ 196 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 588 ]), &(nmpcWorkspace.evGx[ 392 ]), &(nmpcWorkspace.QGx[ 392 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 784 ]), &(nmpcWorkspace.evGx[ 588 ]), &(nmpcWorkspace.QGx[ 588 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 980 ]), &(nmpcWorkspace.evGx[ 784 ]), &(nmpcWorkspace.QGx[ 784 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 1176 ]), &(nmpcWorkspace.evGx[ 980 ]), &(nmpcWorkspace.QGx[ 980 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 1372 ]), &(nmpcWorkspace.evGx[ 1176 ]), &(nmpcWorkspace.QGx[ 1176 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 1568 ]), &(nmpcWorkspace.evGx[ 1372 ]), &(nmpcWorkspace.QGx[ 1372 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 1764 ]), &(nmpcWorkspace.evGx[ 1568 ]), &(nmpcWorkspace.QGx[ 1568 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 1960 ]), &(nmpcWorkspace.evGx[ 1764 ]), &(nmpcWorkspace.QGx[ 1764 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 2156 ]), &(nmpcWorkspace.evGx[ 1960 ]), &(nmpcWorkspace.QGx[ 1960 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 2352 ]), &(nmpcWorkspace.evGx[ 2156 ]), &(nmpcWorkspace.QGx[ 2156 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 2548 ]), &(nmpcWorkspace.evGx[ 2352 ]), &(nmpcWorkspace.QGx[ 2352 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 2744 ]), &(nmpcWorkspace.evGx[ 2548 ]), &(nmpcWorkspace.QGx[ 2548 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 2940 ]), &(nmpcWorkspace.evGx[ 2744 ]), &(nmpcWorkspace.QGx[ 2744 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 3136 ]), &(nmpcWorkspace.evGx[ 2940 ]), &(nmpcWorkspace.QGx[ 2940 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 3332 ]), &(nmpcWorkspace.evGx[ 3136 ]), &(nmpcWorkspace.QGx[ 3136 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 3528 ]), &(nmpcWorkspace.evGx[ 3332 ]), &(nmpcWorkspace.QGx[ 3332 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 3724 ]), &(nmpcWorkspace.evGx[ 3528 ]), &(nmpcWorkspace.QGx[ 3528 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 3920 ]), &(nmpcWorkspace.evGx[ 3724 ]), &(nmpcWorkspace.QGx[ 3724 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 4116 ]), &(nmpcWorkspace.evGx[ 3920 ]), &(nmpcWorkspace.QGx[ 3920 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 4312 ]), &(nmpcWorkspace.evGx[ 4116 ]), &(nmpcWorkspace.QGx[ 4116 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 4508 ]), &(nmpcWorkspace.evGx[ 4312 ]), &(nmpcWorkspace.QGx[ 4312 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 4704 ]), &(nmpcWorkspace.evGx[ 4508 ]), &(nmpcWorkspace.QGx[ 4508 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 4900 ]), &(nmpcWorkspace.evGx[ 4704 ]), &(nmpcWorkspace.QGx[ 4704 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 5096 ]), &(nmpcWorkspace.evGx[ 4900 ]), &(nmpcWorkspace.QGx[ 4900 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 5292 ]), &(nmpcWorkspace.evGx[ 5096 ]), &(nmpcWorkspace.QGx[ 5096 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 5488 ]), &(nmpcWorkspace.evGx[ 5292 ]), &(nmpcWorkspace.QGx[ 5292 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 5684 ]), &(nmpcWorkspace.evGx[ 5488 ]), &(nmpcWorkspace.QGx[ 5488 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 5880 ]), &(nmpcWorkspace.evGx[ 5684 ]), &(nmpcWorkspace.QGx[ 5684 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 6076 ]), &(nmpcWorkspace.evGx[ 5880 ]), &(nmpcWorkspace.QGx[ 5880 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 6272 ]), &(nmpcWorkspace.evGx[ 6076 ]), &(nmpcWorkspace.QGx[ 6076 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 6468 ]), &(nmpcWorkspace.evGx[ 6272 ]), &(nmpcWorkspace.QGx[ 6272 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 6664 ]), &(nmpcWorkspace.evGx[ 6468 ]), &(nmpcWorkspace.QGx[ 6468 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 6860 ]), &(nmpcWorkspace.evGx[ 6664 ]), &(nmpcWorkspace.QGx[ 6664 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 7056 ]), &(nmpcWorkspace.evGx[ 6860 ]), &(nmpcWorkspace.QGx[ 6860 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 7252 ]), &(nmpcWorkspace.evGx[ 7056 ]), &(nmpcWorkspace.QGx[ 7056 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 7448 ]), &(nmpcWorkspace.evGx[ 7252 ]), &(nmpcWorkspace.QGx[ 7252 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 7644 ]), &(nmpcWorkspace.evGx[ 7448 ]), &(nmpcWorkspace.QGx[ 7448 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 7840 ]), &(nmpcWorkspace.evGx[ 7644 ]), &(nmpcWorkspace.QGx[ 7644 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 8036 ]), &(nmpcWorkspace.evGx[ 7840 ]), &(nmpcWorkspace.QGx[ 7840 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 8232 ]), &(nmpcWorkspace.evGx[ 8036 ]), &(nmpcWorkspace.QGx[ 8036 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 8428 ]), &(nmpcWorkspace.evGx[ 8232 ]), &(nmpcWorkspace.QGx[ 8232 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 8624 ]), &(nmpcWorkspace.evGx[ 8428 ]), &(nmpcWorkspace.QGx[ 8428 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 8820 ]), &(nmpcWorkspace.evGx[ 8624 ]), &(nmpcWorkspace.QGx[ 8624 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 9016 ]), &(nmpcWorkspace.evGx[ 8820 ]), &(nmpcWorkspace.QGx[ 8820 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 9212 ]), &(nmpcWorkspace.evGx[ 9016 ]), &(nmpcWorkspace.QGx[ 9016 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 9408 ]), &(nmpcWorkspace.evGx[ 9212 ]), &(nmpcWorkspace.QGx[ 9212 ]) );
nmpc_multGxGx( &(nmpcWorkspace.Q1[ 9604 ]), &(nmpcWorkspace.evGx[ 9408 ]), &(nmpcWorkspace.QGx[ 9408 ]) );
nmpc_multGxGx( nmpcWorkspace.QN1, &(nmpcWorkspace.evGx[ 9604 ]), &(nmpcWorkspace.QGx[ 9604 ]) );

for (lRun1 = 0; lRun1 < 49; ++lRun1)
{
for (lRun2 = 0; lRun2 < lRun1 + 1; ++lRun2)
{
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
nmpc_multGxGu( &(nmpcWorkspace.Q1[ lRun1 * 196 + 196 ]), &(nmpcWorkspace.E[ lRun3 * 56 ]), &(nmpcWorkspace.QE[ lRun3 * 56 ]) );
}
}

for (lRun2 = 0; lRun2 < lRun1 + 1; ++lRun2)
{
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
nmpc_multGxGu( nmpcWorkspace.QN1, &(nmpcWorkspace.E[ lRun3 * 56 ]), &(nmpcWorkspace.QE[ lRun3 * 56 ]) );
}

nmpc_zeroBlockH00(  );
nmpc_multCTQC( nmpcWorkspace.evGx, nmpcWorkspace.QGx );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 196 ]), &(nmpcWorkspace.QGx[ 196 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 392 ]), &(nmpcWorkspace.QGx[ 392 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 588 ]), &(nmpcWorkspace.QGx[ 588 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 784 ]), &(nmpcWorkspace.QGx[ 784 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 980 ]), &(nmpcWorkspace.QGx[ 980 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 1176 ]), &(nmpcWorkspace.QGx[ 1176 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 1372 ]), &(nmpcWorkspace.QGx[ 1372 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 1568 ]), &(nmpcWorkspace.QGx[ 1568 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 1764 ]), &(nmpcWorkspace.QGx[ 1764 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 1960 ]), &(nmpcWorkspace.QGx[ 1960 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 2156 ]), &(nmpcWorkspace.QGx[ 2156 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 2352 ]), &(nmpcWorkspace.QGx[ 2352 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 2548 ]), &(nmpcWorkspace.QGx[ 2548 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 2744 ]), &(nmpcWorkspace.QGx[ 2744 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 2940 ]), &(nmpcWorkspace.QGx[ 2940 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 3136 ]), &(nmpcWorkspace.QGx[ 3136 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 3332 ]), &(nmpcWorkspace.QGx[ 3332 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 3528 ]), &(nmpcWorkspace.QGx[ 3528 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 3724 ]), &(nmpcWorkspace.QGx[ 3724 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 3920 ]), &(nmpcWorkspace.QGx[ 3920 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 4116 ]), &(nmpcWorkspace.QGx[ 4116 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 4312 ]), &(nmpcWorkspace.QGx[ 4312 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 4508 ]), &(nmpcWorkspace.QGx[ 4508 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 4704 ]), &(nmpcWorkspace.QGx[ 4704 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 4900 ]), &(nmpcWorkspace.QGx[ 4900 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 5096 ]), &(nmpcWorkspace.QGx[ 5096 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 5292 ]), &(nmpcWorkspace.QGx[ 5292 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 5488 ]), &(nmpcWorkspace.QGx[ 5488 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 5684 ]), &(nmpcWorkspace.QGx[ 5684 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 5880 ]), &(nmpcWorkspace.QGx[ 5880 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 6076 ]), &(nmpcWorkspace.QGx[ 6076 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 6272 ]), &(nmpcWorkspace.QGx[ 6272 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 6468 ]), &(nmpcWorkspace.QGx[ 6468 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 6664 ]), &(nmpcWorkspace.QGx[ 6664 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 6860 ]), &(nmpcWorkspace.QGx[ 6860 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 7056 ]), &(nmpcWorkspace.QGx[ 7056 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 7252 ]), &(nmpcWorkspace.QGx[ 7252 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 7448 ]), &(nmpcWorkspace.QGx[ 7448 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 7644 ]), &(nmpcWorkspace.QGx[ 7644 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 7840 ]), &(nmpcWorkspace.QGx[ 7840 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 8036 ]), &(nmpcWorkspace.QGx[ 8036 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 8232 ]), &(nmpcWorkspace.QGx[ 8232 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 8428 ]), &(nmpcWorkspace.QGx[ 8428 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 8624 ]), &(nmpcWorkspace.QGx[ 8624 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 8820 ]), &(nmpcWorkspace.QGx[ 8820 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 9016 ]), &(nmpcWorkspace.QGx[ 9016 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 9212 ]), &(nmpcWorkspace.QGx[ 9212 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 9408 ]), &(nmpcWorkspace.QGx[ 9408 ]) );
nmpc_multCTQC( &(nmpcWorkspace.evGx[ 9604 ]), &(nmpcWorkspace.QGx[ 9604 ]) );

for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
nmpc_zeroBlockH10( &(nmpcWorkspace.H10[ lRun1 * 56 ]) );
for (lRun2 = lRun1; lRun2 < 50; ++lRun2)
{
lRun3 = (((lRun2 + 1) * (lRun2)) / (2)) + (lRun1);
nmpc_multQETGx( &(nmpcWorkspace.QE[ lRun3 * 56 ]), &(nmpcWorkspace.evGx[ lRun2 * 196 ]), &(nmpcWorkspace.H10[ lRun1 * 56 ]) );
}
}

for (lRun1 = 0;lRun1 < 14; ++lRun1)
for (lRun2 = 0;lRun2 < 200; ++lRun2)
nmpcWorkspace.H[(lRun1 * 214) + (lRun2 + 14)] = nmpcWorkspace.H10[(lRun2 * 14) + (lRun1)];

for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
nmpc_setBlockH11_R1( lRun1, lRun1, &(nmpcWorkspace.R1[ lRun1 * 16 ]) );
lRun2 = lRun1;
for (lRun3 = lRun1; lRun3 < 50; ++lRun3)
{
lRun4 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun1);
lRun5 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun2);
nmpc_setBlockH11( lRun1, lRun2, &(nmpcWorkspace.E[ lRun4 * 56 ]), &(nmpcWorkspace.QE[ lRun5 * 56 ]) );
}
for (lRun2 = lRun1 + 1; lRun2 < 50; ++lRun2)
{
nmpc_zeroBlockH11( lRun1, lRun2 );
for (lRun3 = lRun2; lRun3 < 50; ++lRun3)
{
lRun4 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun1);
lRun5 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun2);
nmpc_setBlockH11( lRun1, lRun2, &(nmpcWorkspace.E[ lRun4 * 56 ]), &(nmpcWorkspace.QE[ lRun5 * 56 ]) );
}
}
}

for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
for (lRun2 = 0; lRun2 < lRun1; ++lRun2)
{
nmpc_copyHTH( lRun1, lRun2 );
}
}

for (lRun1 = 0;lRun1 < 200; ++lRun1)
for (lRun2 = 0;lRun2 < 14; ++lRun2)
nmpcWorkspace.H[(lRun1 * 214 + 2996) + (lRun2)] = nmpcWorkspace.H10[(lRun1 * 14) + (lRun2)];

nmpc_multQ1d( &(nmpcWorkspace.Q1[ 196 ]), nmpcWorkspace.d, nmpcWorkspace.Qd );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 392 ]), &(nmpcWorkspace.d[ 14 ]), &(nmpcWorkspace.Qd[ 14 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 588 ]), &(nmpcWorkspace.d[ 28 ]), &(nmpcWorkspace.Qd[ 28 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 784 ]), &(nmpcWorkspace.d[ 42 ]), &(nmpcWorkspace.Qd[ 42 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 980 ]), &(nmpcWorkspace.d[ 56 ]), &(nmpcWorkspace.Qd[ 56 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 1176 ]), &(nmpcWorkspace.d[ 70 ]), &(nmpcWorkspace.Qd[ 70 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 1372 ]), &(nmpcWorkspace.d[ 84 ]), &(nmpcWorkspace.Qd[ 84 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 1568 ]), &(nmpcWorkspace.d[ 98 ]), &(nmpcWorkspace.Qd[ 98 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 1764 ]), &(nmpcWorkspace.d[ 112 ]), &(nmpcWorkspace.Qd[ 112 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 1960 ]), &(nmpcWorkspace.d[ 126 ]), &(nmpcWorkspace.Qd[ 126 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 2156 ]), &(nmpcWorkspace.d[ 140 ]), &(nmpcWorkspace.Qd[ 140 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 2352 ]), &(nmpcWorkspace.d[ 154 ]), &(nmpcWorkspace.Qd[ 154 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 2548 ]), &(nmpcWorkspace.d[ 168 ]), &(nmpcWorkspace.Qd[ 168 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 2744 ]), &(nmpcWorkspace.d[ 182 ]), &(nmpcWorkspace.Qd[ 182 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 2940 ]), &(nmpcWorkspace.d[ 196 ]), &(nmpcWorkspace.Qd[ 196 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 3136 ]), &(nmpcWorkspace.d[ 210 ]), &(nmpcWorkspace.Qd[ 210 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 3332 ]), &(nmpcWorkspace.d[ 224 ]), &(nmpcWorkspace.Qd[ 224 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 3528 ]), &(nmpcWorkspace.d[ 238 ]), &(nmpcWorkspace.Qd[ 238 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 3724 ]), &(nmpcWorkspace.d[ 252 ]), &(nmpcWorkspace.Qd[ 252 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 3920 ]), &(nmpcWorkspace.d[ 266 ]), &(nmpcWorkspace.Qd[ 266 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 4116 ]), &(nmpcWorkspace.d[ 280 ]), &(nmpcWorkspace.Qd[ 280 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 4312 ]), &(nmpcWorkspace.d[ 294 ]), &(nmpcWorkspace.Qd[ 294 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 4508 ]), &(nmpcWorkspace.d[ 308 ]), &(nmpcWorkspace.Qd[ 308 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 4704 ]), &(nmpcWorkspace.d[ 322 ]), &(nmpcWorkspace.Qd[ 322 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 4900 ]), &(nmpcWorkspace.d[ 336 ]), &(nmpcWorkspace.Qd[ 336 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 5096 ]), &(nmpcWorkspace.d[ 350 ]), &(nmpcWorkspace.Qd[ 350 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 5292 ]), &(nmpcWorkspace.d[ 364 ]), &(nmpcWorkspace.Qd[ 364 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 5488 ]), &(nmpcWorkspace.d[ 378 ]), &(nmpcWorkspace.Qd[ 378 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 5684 ]), &(nmpcWorkspace.d[ 392 ]), &(nmpcWorkspace.Qd[ 392 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 5880 ]), &(nmpcWorkspace.d[ 406 ]), &(nmpcWorkspace.Qd[ 406 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 6076 ]), &(nmpcWorkspace.d[ 420 ]), &(nmpcWorkspace.Qd[ 420 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 6272 ]), &(nmpcWorkspace.d[ 434 ]), &(nmpcWorkspace.Qd[ 434 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 6468 ]), &(nmpcWorkspace.d[ 448 ]), &(nmpcWorkspace.Qd[ 448 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 6664 ]), &(nmpcWorkspace.d[ 462 ]), &(nmpcWorkspace.Qd[ 462 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 6860 ]), &(nmpcWorkspace.d[ 476 ]), &(nmpcWorkspace.Qd[ 476 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 7056 ]), &(nmpcWorkspace.d[ 490 ]), &(nmpcWorkspace.Qd[ 490 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 7252 ]), &(nmpcWorkspace.d[ 504 ]), &(nmpcWorkspace.Qd[ 504 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 7448 ]), &(nmpcWorkspace.d[ 518 ]), &(nmpcWorkspace.Qd[ 518 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 7644 ]), &(nmpcWorkspace.d[ 532 ]), &(nmpcWorkspace.Qd[ 532 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 7840 ]), &(nmpcWorkspace.d[ 546 ]), &(nmpcWorkspace.Qd[ 546 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 8036 ]), &(nmpcWorkspace.d[ 560 ]), &(nmpcWorkspace.Qd[ 560 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 8232 ]), &(nmpcWorkspace.d[ 574 ]), &(nmpcWorkspace.Qd[ 574 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 8428 ]), &(nmpcWorkspace.d[ 588 ]), &(nmpcWorkspace.Qd[ 588 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 8624 ]), &(nmpcWorkspace.d[ 602 ]), &(nmpcWorkspace.Qd[ 602 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 8820 ]), &(nmpcWorkspace.d[ 616 ]), &(nmpcWorkspace.Qd[ 616 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 9016 ]), &(nmpcWorkspace.d[ 630 ]), &(nmpcWorkspace.Qd[ 630 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 9212 ]), &(nmpcWorkspace.d[ 644 ]), &(nmpcWorkspace.Qd[ 644 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 9408 ]), &(nmpcWorkspace.d[ 658 ]), &(nmpcWorkspace.Qd[ 658 ]) );
nmpc_multQ1d( &(nmpcWorkspace.Q1[ 9604 ]), &(nmpcWorkspace.d[ 672 ]), &(nmpcWorkspace.Qd[ 672 ]) );
nmpc_multQN1d( nmpcWorkspace.QN1, &(nmpcWorkspace.d[ 686 ]), &(nmpcWorkspace.Qd[ 686 ]) );

nmpc_macCTSlx( nmpcWorkspace.evGx, nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 196 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 392 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 588 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 784 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 980 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 1176 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 1372 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 1568 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 1764 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 1960 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 2156 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 2352 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 2548 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 2744 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 2940 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 3136 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 3332 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 3528 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 3724 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 3920 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 4116 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 4312 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 4508 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 4704 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 4900 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 5096 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 5292 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 5488 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 5684 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 5880 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 6076 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 6272 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 6468 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 6664 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 6860 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 7056 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 7252 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 7448 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 7644 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 7840 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 8036 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 8232 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 8428 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 8624 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 8820 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 9016 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 9212 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 9408 ]), nmpcWorkspace.g );
nmpc_macCTSlx( &(nmpcWorkspace.evGx[ 9604 ]), nmpcWorkspace.g );
for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
for (lRun2 = lRun1; lRun2 < 50; ++lRun2)
{
lRun3 = (((lRun2 + 1) * (lRun2)) / (2)) + (lRun1);
nmpc_macETSlu( &(nmpcWorkspace.QE[ lRun3 * 56 ]), &(nmpcWorkspace.g[ lRun1 * 4 + 14 ]) );
}
}
nmpcWorkspace.lb[14] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[0];
nmpcWorkspace.lb[15] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[1];
nmpcWorkspace.lb[16] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[2];
nmpcWorkspace.lb[17] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[3];
nmpcWorkspace.lb[18] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[4];
nmpcWorkspace.lb[19] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[5];
nmpcWorkspace.lb[20] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[6];
nmpcWorkspace.lb[21] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[7];
nmpcWorkspace.lb[22] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[8];
nmpcWorkspace.lb[23] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[9];
nmpcWorkspace.lb[24] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[10];
nmpcWorkspace.lb[25] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[11];
nmpcWorkspace.lb[26] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[12];
nmpcWorkspace.lb[27] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[13];
nmpcWorkspace.lb[28] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[14];
nmpcWorkspace.lb[29] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[15];
nmpcWorkspace.lb[30] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[16];
nmpcWorkspace.lb[31] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[17];
nmpcWorkspace.lb[32] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[18];
nmpcWorkspace.lb[33] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[19];
nmpcWorkspace.lb[34] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[20];
nmpcWorkspace.lb[35] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[21];
nmpcWorkspace.lb[36] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[22];
nmpcWorkspace.lb[37] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[23];
nmpcWorkspace.lb[38] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[24];
nmpcWorkspace.lb[39] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[25];
nmpcWorkspace.lb[40] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[26];
nmpcWorkspace.lb[41] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[27];
nmpcWorkspace.lb[42] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[28];
nmpcWorkspace.lb[43] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[29];
nmpcWorkspace.lb[44] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[30];
nmpcWorkspace.lb[45] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[31];
nmpcWorkspace.lb[46] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[32];
nmpcWorkspace.lb[47] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[33];
nmpcWorkspace.lb[48] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[34];
nmpcWorkspace.lb[49] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[35];
nmpcWorkspace.lb[50] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[36];
nmpcWorkspace.lb[51] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[37];
nmpcWorkspace.lb[52] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[38];
nmpcWorkspace.lb[53] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[39];
nmpcWorkspace.lb[54] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[40];
nmpcWorkspace.lb[55] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[41];
nmpcWorkspace.lb[56] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[42];
nmpcWorkspace.lb[57] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[43];
nmpcWorkspace.lb[58] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[44];
nmpcWorkspace.lb[59] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[45];
nmpcWorkspace.lb[60] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[46];
nmpcWorkspace.lb[61] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[47];
nmpcWorkspace.lb[62] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[48];
nmpcWorkspace.lb[63] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[49];
nmpcWorkspace.lb[64] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[50];
nmpcWorkspace.lb[65] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[51];
nmpcWorkspace.lb[66] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[52];
nmpcWorkspace.lb[67] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[53];
nmpcWorkspace.lb[68] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[54];
nmpcWorkspace.lb[69] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[55];
nmpcWorkspace.lb[70] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[56];
nmpcWorkspace.lb[71] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[57];
nmpcWorkspace.lb[72] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[58];
nmpcWorkspace.lb[73] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[59];
nmpcWorkspace.lb[74] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[60];
nmpcWorkspace.lb[75] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[61];
nmpcWorkspace.lb[76] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[62];
nmpcWorkspace.lb[77] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[63];
nmpcWorkspace.lb[78] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[64];
nmpcWorkspace.lb[79] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[65];
nmpcWorkspace.lb[80] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[66];
nmpcWorkspace.lb[81] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[67];
nmpcWorkspace.lb[82] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[68];
nmpcWorkspace.lb[83] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[69];
nmpcWorkspace.lb[84] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[70];
nmpcWorkspace.lb[85] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[71];
nmpcWorkspace.lb[86] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[72];
nmpcWorkspace.lb[87] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[73];
nmpcWorkspace.lb[88] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[74];
nmpcWorkspace.lb[89] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[75];
nmpcWorkspace.lb[90] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[76];
nmpcWorkspace.lb[91] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[77];
nmpcWorkspace.lb[92] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[78];
nmpcWorkspace.lb[93] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[79];
nmpcWorkspace.lb[94] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[80];
nmpcWorkspace.lb[95] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[81];
nmpcWorkspace.lb[96] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[82];
nmpcWorkspace.lb[97] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[83];
nmpcWorkspace.lb[98] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[84];
nmpcWorkspace.lb[99] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[85];
nmpcWorkspace.lb[100] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[86];
nmpcWorkspace.lb[101] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[87];
nmpcWorkspace.lb[102] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[88];
nmpcWorkspace.lb[103] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[89];
nmpcWorkspace.lb[104] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[90];
nmpcWorkspace.lb[105] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[91];
nmpcWorkspace.lb[106] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[92];
nmpcWorkspace.lb[107] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[93];
nmpcWorkspace.lb[108] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[94];
nmpcWorkspace.lb[109] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[95];
nmpcWorkspace.lb[110] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[96];
nmpcWorkspace.lb[111] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[97];
nmpcWorkspace.lb[112] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[98];
nmpcWorkspace.lb[113] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[99];
nmpcWorkspace.lb[114] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[100];
nmpcWorkspace.lb[115] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[101];
nmpcWorkspace.lb[116] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[102];
nmpcWorkspace.lb[117] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[103];
nmpcWorkspace.lb[118] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[104];
nmpcWorkspace.lb[119] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[105];
nmpcWorkspace.lb[120] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[106];
nmpcWorkspace.lb[121] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[107];
nmpcWorkspace.lb[122] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[108];
nmpcWorkspace.lb[123] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[109];
nmpcWorkspace.lb[124] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[110];
nmpcWorkspace.lb[125] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[111];
nmpcWorkspace.lb[126] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[112];
nmpcWorkspace.lb[127] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[113];
nmpcWorkspace.lb[128] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[114];
nmpcWorkspace.lb[129] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[115];
nmpcWorkspace.lb[130] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[116];
nmpcWorkspace.lb[131] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[117];
nmpcWorkspace.lb[132] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[118];
nmpcWorkspace.lb[133] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[119];
nmpcWorkspace.lb[134] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[120];
nmpcWorkspace.lb[135] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[121];
nmpcWorkspace.lb[136] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[122];
nmpcWorkspace.lb[137] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[123];
nmpcWorkspace.lb[138] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[124];
nmpcWorkspace.lb[139] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[125];
nmpcWorkspace.lb[140] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[126];
nmpcWorkspace.lb[141] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[127];
nmpcWorkspace.lb[142] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[128];
nmpcWorkspace.lb[143] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[129];
nmpcWorkspace.lb[144] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[130];
nmpcWorkspace.lb[145] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[131];
nmpcWorkspace.lb[146] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[132];
nmpcWorkspace.lb[147] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[133];
nmpcWorkspace.lb[148] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[134];
nmpcWorkspace.lb[149] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[135];
nmpcWorkspace.lb[150] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[136];
nmpcWorkspace.lb[151] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[137];
nmpcWorkspace.lb[152] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[138];
nmpcWorkspace.lb[153] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[139];
nmpcWorkspace.lb[154] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[140];
nmpcWorkspace.lb[155] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[141];
nmpcWorkspace.lb[156] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[142];
nmpcWorkspace.lb[157] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[143];
nmpcWorkspace.lb[158] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[144];
nmpcWorkspace.lb[159] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[145];
nmpcWorkspace.lb[160] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[146];
nmpcWorkspace.lb[161] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[147];
nmpcWorkspace.lb[162] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[148];
nmpcWorkspace.lb[163] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[149];
nmpcWorkspace.lb[164] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[150];
nmpcWorkspace.lb[165] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[151];
nmpcWorkspace.lb[166] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[152];
nmpcWorkspace.lb[167] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[153];
nmpcWorkspace.lb[168] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[154];
nmpcWorkspace.lb[169] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[155];
nmpcWorkspace.lb[170] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[156];
nmpcWorkspace.lb[171] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[157];
nmpcWorkspace.lb[172] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[158];
nmpcWorkspace.lb[173] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[159];
nmpcWorkspace.lb[174] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[160];
nmpcWorkspace.lb[175] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[161];
nmpcWorkspace.lb[176] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[162];
nmpcWorkspace.lb[177] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[163];
nmpcWorkspace.lb[178] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[164];
nmpcWorkspace.lb[179] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[165];
nmpcWorkspace.lb[180] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[166];
nmpcWorkspace.lb[181] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[167];
nmpcWorkspace.lb[182] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[168];
nmpcWorkspace.lb[183] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[169];
nmpcWorkspace.lb[184] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[170];
nmpcWorkspace.lb[185] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[171];
nmpcWorkspace.lb[186] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[172];
nmpcWorkspace.lb[187] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[173];
nmpcWorkspace.lb[188] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[174];
nmpcWorkspace.lb[189] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[175];
nmpcWorkspace.lb[190] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[176];
nmpcWorkspace.lb[191] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[177];
nmpcWorkspace.lb[192] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[178];
nmpcWorkspace.lb[193] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[179];
nmpcWorkspace.lb[194] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[180];
nmpcWorkspace.lb[195] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[181];
nmpcWorkspace.lb[196] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[182];
nmpcWorkspace.lb[197] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[183];
nmpcWorkspace.lb[198] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[184];
nmpcWorkspace.lb[199] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[185];
nmpcWorkspace.lb[200] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[186];
nmpcWorkspace.lb[201] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[187];
nmpcWorkspace.lb[202] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[188];
nmpcWorkspace.lb[203] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[189];
nmpcWorkspace.lb[204] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[190];
nmpcWorkspace.lb[205] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[191];
nmpcWorkspace.lb[206] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[192];
nmpcWorkspace.lb[207] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[193];
nmpcWorkspace.lb[208] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[194];
nmpcWorkspace.lb[209] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[195];
nmpcWorkspace.lb[210] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[196];
nmpcWorkspace.lb[211] = (real_t)-8.0000000000000000e+01 - nmpcVariables.u[197];
nmpcWorkspace.lb[212] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[198];
nmpcWorkspace.lb[213] = (real_t)-1.6000000000000000e+02 - nmpcVariables.u[199];
nmpcWorkspace.ub[14] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[0];
nmpcWorkspace.ub[15] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[1];
nmpcWorkspace.ub[16] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[2];
nmpcWorkspace.ub[17] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[3];
nmpcWorkspace.ub[18] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[4];
nmpcWorkspace.ub[19] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[5];
nmpcWorkspace.ub[20] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[6];
nmpcWorkspace.ub[21] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[7];
nmpcWorkspace.ub[22] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[8];
nmpcWorkspace.ub[23] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[9];
nmpcWorkspace.ub[24] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[10];
nmpcWorkspace.ub[25] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[11];
nmpcWorkspace.ub[26] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[12];
nmpcWorkspace.ub[27] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[13];
nmpcWorkspace.ub[28] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[14];
nmpcWorkspace.ub[29] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[15];
nmpcWorkspace.ub[30] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[16];
nmpcWorkspace.ub[31] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[17];
nmpcWorkspace.ub[32] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[18];
nmpcWorkspace.ub[33] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[19];
nmpcWorkspace.ub[34] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[20];
nmpcWorkspace.ub[35] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[21];
nmpcWorkspace.ub[36] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[22];
nmpcWorkspace.ub[37] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[23];
nmpcWorkspace.ub[38] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[24];
nmpcWorkspace.ub[39] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[25];
nmpcWorkspace.ub[40] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[26];
nmpcWorkspace.ub[41] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[27];
nmpcWorkspace.ub[42] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[28];
nmpcWorkspace.ub[43] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[29];
nmpcWorkspace.ub[44] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[30];
nmpcWorkspace.ub[45] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[31];
nmpcWorkspace.ub[46] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[32];
nmpcWorkspace.ub[47] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[33];
nmpcWorkspace.ub[48] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[34];
nmpcWorkspace.ub[49] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[35];
nmpcWorkspace.ub[50] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[36];
nmpcWorkspace.ub[51] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[37];
nmpcWorkspace.ub[52] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[38];
nmpcWorkspace.ub[53] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[39];
nmpcWorkspace.ub[54] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[40];
nmpcWorkspace.ub[55] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[41];
nmpcWorkspace.ub[56] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[42];
nmpcWorkspace.ub[57] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[43];
nmpcWorkspace.ub[58] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[44];
nmpcWorkspace.ub[59] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[45];
nmpcWorkspace.ub[60] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[46];
nmpcWorkspace.ub[61] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[47];
nmpcWorkspace.ub[62] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[48];
nmpcWorkspace.ub[63] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[49];
nmpcWorkspace.ub[64] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[50];
nmpcWorkspace.ub[65] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[51];
nmpcWorkspace.ub[66] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[52];
nmpcWorkspace.ub[67] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[53];
nmpcWorkspace.ub[68] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[54];
nmpcWorkspace.ub[69] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[55];
nmpcWorkspace.ub[70] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[56];
nmpcWorkspace.ub[71] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[57];
nmpcWorkspace.ub[72] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[58];
nmpcWorkspace.ub[73] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[59];
nmpcWorkspace.ub[74] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[60];
nmpcWorkspace.ub[75] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[61];
nmpcWorkspace.ub[76] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[62];
nmpcWorkspace.ub[77] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[63];
nmpcWorkspace.ub[78] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[64];
nmpcWorkspace.ub[79] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[65];
nmpcWorkspace.ub[80] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[66];
nmpcWorkspace.ub[81] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[67];
nmpcWorkspace.ub[82] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[68];
nmpcWorkspace.ub[83] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[69];
nmpcWorkspace.ub[84] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[70];
nmpcWorkspace.ub[85] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[71];
nmpcWorkspace.ub[86] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[72];
nmpcWorkspace.ub[87] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[73];
nmpcWorkspace.ub[88] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[74];
nmpcWorkspace.ub[89] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[75];
nmpcWorkspace.ub[90] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[76];
nmpcWorkspace.ub[91] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[77];
nmpcWorkspace.ub[92] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[78];
nmpcWorkspace.ub[93] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[79];
nmpcWorkspace.ub[94] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[80];
nmpcWorkspace.ub[95] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[81];
nmpcWorkspace.ub[96] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[82];
nmpcWorkspace.ub[97] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[83];
nmpcWorkspace.ub[98] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[84];
nmpcWorkspace.ub[99] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[85];
nmpcWorkspace.ub[100] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[86];
nmpcWorkspace.ub[101] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[87];
nmpcWorkspace.ub[102] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[88];
nmpcWorkspace.ub[103] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[89];
nmpcWorkspace.ub[104] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[90];
nmpcWorkspace.ub[105] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[91];
nmpcWorkspace.ub[106] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[92];
nmpcWorkspace.ub[107] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[93];
nmpcWorkspace.ub[108] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[94];
nmpcWorkspace.ub[109] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[95];
nmpcWorkspace.ub[110] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[96];
nmpcWorkspace.ub[111] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[97];
nmpcWorkspace.ub[112] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[98];
nmpcWorkspace.ub[113] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[99];
nmpcWorkspace.ub[114] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[100];
nmpcWorkspace.ub[115] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[101];
nmpcWorkspace.ub[116] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[102];
nmpcWorkspace.ub[117] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[103];
nmpcWorkspace.ub[118] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[104];
nmpcWorkspace.ub[119] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[105];
nmpcWorkspace.ub[120] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[106];
nmpcWorkspace.ub[121] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[107];
nmpcWorkspace.ub[122] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[108];
nmpcWorkspace.ub[123] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[109];
nmpcWorkspace.ub[124] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[110];
nmpcWorkspace.ub[125] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[111];
nmpcWorkspace.ub[126] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[112];
nmpcWorkspace.ub[127] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[113];
nmpcWorkspace.ub[128] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[114];
nmpcWorkspace.ub[129] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[115];
nmpcWorkspace.ub[130] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[116];
nmpcWorkspace.ub[131] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[117];
nmpcWorkspace.ub[132] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[118];
nmpcWorkspace.ub[133] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[119];
nmpcWorkspace.ub[134] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[120];
nmpcWorkspace.ub[135] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[121];
nmpcWorkspace.ub[136] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[122];
nmpcWorkspace.ub[137] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[123];
nmpcWorkspace.ub[138] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[124];
nmpcWorkspace.ub[139] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[125];
nmpcWorkspace.ub[140] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[126];
nmpcWorkspace.ub[141] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[127];
nmpcWorkspace.ub[142] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[128];
nmpcWorkspace.ub[143] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[129];
nmpcWorkspace.ub[144] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[130];
nmpcWorkspace.ub[145] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[131];
nmpcWorkspace.ub[146] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[132];
nmpcWorkspace.ub[147] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[133];
nmpcWorkspace.ub[148] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[134];
nmpcWorkspace.ub[149] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[135];
nmpcWorkspace.ub[150] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[136];
nmpcWorkspace.ub[151] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[137];
nmpcWorkspace.ub[152] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[138];
nmpcWorkspace.ub[153] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[139];
nmpcWorkspace.ub[154] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[140];
nmpcWorkspace.ub[155] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[141];
nmpcWorkspace.ub[156] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[142];
nmpcWorkspace.ub[157] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[143];
nmpcWorkspace.ub[158] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[144];
nmpcWorkspace.ub[159] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[145];
nmpcWorkspace.ub[160] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[146];
nmpcWorkspace.ub[161] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[147];
nmpcWorkspace.ub[162] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[148];
nmpcWorkspace.ub[163] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[149];
nmpcWorkspace.ub[164] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[150];
nmpcWorkspace.ub[165] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[151];
nmpcWorkspace.ub[166] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[152];
nmpcWorkspace.ub[167] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[153];
nmpcWorkspace.ub[168] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[154];
nmpcWorkspace.ub[169] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[155];
nmpcWorkspace.ub[170] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[156];
nmpcWorkspace.ub[171] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[157];
nmpcWorkspace.ub[172] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[158];
nmpcWorkspace.ub[173] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[159];
nmpcWorkspace.ub[174] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[160];
nmpcWorkspace.ub[175] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[161];
nmpcWorkspace.ub[176] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[162];
nmpcWorkspace.ub[177] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[163];
nmpcWorkspace.ub[178] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[164];
nmpcWorkspace.ub[179] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[165];
nmpcWorkspace.ub[180] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[166];
nmpcWorkspace.ub[181] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[167];
nmpcWorkspace.ub[182] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[168];
nmpcWorkspace.ub[183] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[169];
nmpcWorkspace.ub[184] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[170];
nmpcWorkspace.ub[185] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[171];
nmpcWorkspace.ub[186] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[172];
nmpcWorkspace.ub[187] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[173];
nmpcWorkspace.ub[188] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[174];
nmpcWorkspace.ub[189] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[175];
nmpcWorkspace.ub[190] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[176];
nmpcWorkspace.ub[191] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[177];
nmpcWorkspace.ub[192] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[178];
nmpcWorkspace.ub[193] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[179];
nmpcWorkspace.ub[194] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[180];
nmpcWorkspace.ub[195] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[181];
nmpcWorkspace.ub[196] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[182];
nmpcWorkspace.ub[197] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[183];
nmpcWorkspace.ub[198] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[184];
nmpcWorkspace.ub[199] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[185];
nmpcWorkspace.ub[200] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[186];
nmpcWorkspace.ub[201] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[187];
nmpcWorkspace.ub[202] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[188];
nmpcWorkspace.ub[203] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[189];
nmpcWorkspace.ub[204] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[190];
nmpcWorkspace.ub[205] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[191];
nmpcWorkspace.ub[206] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[192];
nmpcWorkspace.ub[207] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[193];
nmpcWorkspace.ub[208] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[194];
nmpcWorkspace.ub[209] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[195];
nmpcWorkspace.ub[210] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[196];
nmpcWorkspace.ub[211] = (real_t)8.0000000000000000e+01 - nmpcVariables.u[197];
nmpcWorkspace.ub[212] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[198];
nmpcWorkspace.ub[213] = (real_t)1.6000000000000000e+02 - nmpcVariables.u[199];

}

void nmpc_condenseFdb(  )
{
int lRun1;
int lRun2;
int lRun3;
int lRun4;
int lRun5;
nmpcWorkspace.Dx0[0] = nmpcVariables.x0[0] - nmpcVariables.x[0];
nmpcWorkspace.Dx0[1] = nmpcVariables.x0[1] - nmpcVariables.x[1];
nmpcWorkspace.Dx0[2] = nmpcVariables.x0[2] - nmpcVariables.x[2];
nmpcWorkspace.Dx0[3] = nmpcVariables.x0[3] - nmpcVariables.x[3];
nmpcWorkspace.Dx0[4] = nmpcVariables.x0[4] - nmpcVariables.x[4];
nmpcWorkspace.Dx0[5] = nmpcVariables.x0[5] - nmpcVariables.x[5];
nmpcWorkspace.Dx0[6] = nmpcVariables.x0[6] - nmpcVariables.x[6];
nmpcWorkspace.Dx0[7] = nmpcVariables.x0[7] - nmpcVariables.x[7];
nmpcWorkspace.Dx0[8] = nmpcVariables.x0[8] - nmpcVariables.x[8];
nmpcWorkspace.Dx0[9] = nmpcVariables.x0[9] - nmpcVariables.x[9];
nmpcWorkspace.Dx0[10] = nmpcVariables.x0[10] - nmpcVariables.x[10];
nmpcWorkspace.Dx0[11] = nmpcVariables.x0[11] - nmpcVariables.x[11];
nmpcWorkspace.Dx0[12] = nmpcVariables.x0[12] - nmpcVariables.x[12];
nmpcWorkspace.Dx0[13] = nmpcVariables.x0[13] - nmpcVariables.x[13];

for (lRun2 = 0; lRun2 < 650; ++lRun2)
nmpcWorkspace.Dy[lRun2] -= nmpcVariables.y[lRun2];

nmpcWorkspace.DyN[0] -= nmpcVariables.yN[0];
nmpcWorkspace.DyN[1] -= nmpcVariables.yN[1];
nmpcWorkspace.DyN[2] -= nmpcVariables.yN[2];
nmpcWorkspace.DyN[3] -= nmpcVariables.yN[3];
nmpcWorkspace.DyN[4] -= nmpcVariables.yN[4];
nmpcWorkspace.DyN[5] -= nmpcVariables.yN[5];
nmpcWorkspace.DyN[6] -= nmpcVariables.yN[6];
nmpcWorkspace.DyN[7] -= nmpcVariables.yN[7];

nmpc_multRDy( nmpcWorkspace.R2, nmpcWorkspace.Dy, &(nmpcWorkspace.g[ 14 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 52 ]), &(nmpcWorkspace.Dy[ 13 ]), &(nmpcWorkspace.g[ 18 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 104 ]), &(nmpcWorkspace.Dy[ 26 ]), &(nmpcWorkspace.g[ 22 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 156 ]), &(nmpcWorkspace.Dy[ 39 ]), &(nmpcWorkspace.g[ 26 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 208 ]), &(nmpcWorkspace.Dy[ 52 ]), &(nmpcWorkspace.g[ 30 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 260 ]), &(nmpcWorkspace.Dy[ 65 ]), &(nmpcWorkspace.g[ 34 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 312 ]), &(nmpcWorkspace.Dy[ 78 ]), &(nmpcWorkspace.g[ 38 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 364 ]), &(nmpcWorkspace.Dy[ 91 ]), &(nmpcWorkspace.g[ 42 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 416 ]), &(nmpcWorkspace.Dy[ 104 ]), &(nmpcWorkspace.g[ 46 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 468 ]), &(nmpcWorkspace.Dy[ 117 ]), &(nmpcWorkspace.g[ 50 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 520 ]), &(nmpcWorkspace.Dy[ 130 ]), &(nmpcWorkspace.g[ 54 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 572 ]), &(nmpcWorkspace.Dy[ 143 ]), &(nmpcWorkspace.g[ 58 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 624 ]), &(nmpcWorkspace.Dy[ 156 ]), &(nmpcWorkspace.g[ 62 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 676 ]), &(nmpcWorkspace.Dy[ 169 ]), &(nmpcWorkspace.g[ 66 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 728 ]), &(nmpcWorkspace.Dy[ 182 ]), &(nmpcWorkspace.g[ 70 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 780 ]), &(nmpcWorkspace.Dy[ 195 ]), &(nmpcWorkspace.g[ 74 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 832 ]), &(nmpcWorkspace.Dy[ 208 ]), &(nmpcWorkspace.g[ 78 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 884 ]), &(nmpcWorkspace.Dy[ 221 ]), &(nmpcWorkspace.g[ 82 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 936 ]), &(nmpcWorkspace.Dy[ 234 ]), &(nmpcWorkspace.g[ 86 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 988 ]), &(nmpcWorkspace.Dy[ 247 ]), &(nmpcWorkspace.g[ 90 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1040 ]), &(nmpcWorkspace.Dy[ 260 ]), &(nmpcWorkspace.g[ 94 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1092 ]), &(nmpcWorkspace.Dy[ 273 ]), &(nmpcWorkspace.g[ 98 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1144 ]), &(nmpcWorkspace.Dy[ 286 ]), &(nmpcWorkspace.g[ 102 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1196 ]), &(nmpcWorkspace.Dy[ 299 ]), &(nmpcWorkspace.g[ 106 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1248 ]), &(nmpcWorkspace.Dy[ 312 ]), &(nmpcWorkspace.g[ 110 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1300 ]), &(nmpcWorkspace.Dy[ 325 ]), &(nmpcWorkspace.g[ 114 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1352 ]), &(nmpcWorkspace.Dy[ 338 ]), &(nmpcWorkspace.g[ 118 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1404 ]), &(nmpcWorkspace.Dy[ 351 ]), &(nmpcWorkspace.g[ 122 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1456 ]), &(nmpcWorkspace.Dy[ 364 ]), &(nmpcWorkspace.g[ 126 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1508 ]), &(nmpcWorkspace.Dy[ 377 ]), &(nmpcWorkspace.g[ 130 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1560 ]), &(nmpcWorkspace.Dy[ 390 ]), &(nmpcWorkspace.g[ 134 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1612 ]), &(nmpcWorkspace.Dy[ 403 ]), &(nmpcWorkspace.g[ 138 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1664 ]), &(nmpcWorkspace.Dy[ 416 ]), &(nmpcWorkspace.g[ 142 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1716 ]), &(nmpcWorkspace.Dy[ 429 ]), &(nmpcWorkspace.g[ 146 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1768 ]), &(nmpcWorkspace.Dy[ 442 ]), &(nmpcWorkspace.g[ 150 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1820 ]), &(nmpcWorkspace.Dy[ 455 ]), &(nmpcWorkspace.g[ 154 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1872 ]), &(nmpcWorkspace.Dy[ 468 ]), &(nmpcWorkspace.g[ 158 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1924 ]), &(nmpcWorkspace.Dy[ 481 ]), &(nmpcWorkspace.g[ 162 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 1976 ]), &(nmpcWorkspace.Dy[ 494 ]), &(nmpcWorkspace.g[ 166 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2028 ]), &(nmpcWorkspace.Dy[ 507 ]), &(nmpcWorkspace.g[ 170 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2080 ]), &(nmpcWorkspace.Dy[ 520 ]), &(nmpcWorkspace.g[ 174 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2132 ]), &(nmpcWorkspace.Dy[ 533 ]), &(nmpcWorkspace.g[ 178 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2184 ]), &(nmpcWorkspace.Dy[ 546 ]), &(nmpcWorkspace.g[ 182 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2236 ]), &(nmpcWorkspace.Dy[ 559 ]), &(nmpcWorkspace.g[ 186 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2288 ]), &(nmpcWorkspace.Dy[ 572 ]), &(nmpcWorkspace.g[ 190 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2340 ]), &(nmpcWorkspace.Dy[ 585 ]), &(nmpcWorkspace.g[ 194 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2392 ]), &(nmpcWorkspace.Dy[ 598 ]), &(nmpcWorkspace.g[ 198 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2444 ]), &(nmpcWorkspace.Dy[ 611 ]), &(nmpcWorkspace.g[ 202 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2496 ]), &(nmpcWorkspace.Dy[ 624 ]), &(nmpcWorkspace.g[ 206 ]) );
nmpc_multRDy( &(nmpcWorkspace.R2[ 2548 ]), &(nmpcWorkspace.Dy[ 637 ]), &(nmpcWorkspace.g[ 210 ]) );

nmpc_multQDy( nmpcWorkspace.Q2, nmpcWorkspace.Dy, nmpcWorkspace.QDy );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 182 ]), &(nmpcWorkspace.Dy[ 13 ]), &(nmpcWorkspace.QDy[ 14 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 364 ]), &(nmpcWorkspace.Dy[ 26 ]), &(nmpcWorkspace.QDy[ 28 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 546 ]), &(nmpcWorkspace.Dy[ 39 ]), &(nmpcWorkspace.QDy[ 42 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 728 ]), &(nmpcWorkspace.Dy[ 52 ]), &(nmpcWorkspace.QDy[ 56 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 910 ]), &(nmpcWorkspace.Dy[ 65 ]), &(nmpcWorkspace.QDy[ 70 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 1092 ]), &(nmpcWorkspace.Dy[ 78 ]), &(nmpcWorkspace.QDy[ 84 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 1274 ]), &(nmpcWorkspace.Dy[ 91 ]), &(nmpcWorkspace.QDy[ 98 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 1456 ]), &(nmpcWorkspace.Dy[ 104 ]), &(nmpcWorkspace.QDy[ 112 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 1638 ]), &(nmpcWorkspace.Dy[ 117 ]), &(nmpcWorkspace.QDy[ 126 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 1820 ]), &(nmpcWorkspace.Dy[ 130 ]), &(nmpcWorkspace.QDy[ 140 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 2002 ]), &(nmpcWorkspace.Dy[ 143 ]), &(nmpcWorkspace.QDy[ 154 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 2184 ]), &(nmpcWorkspace.Dy[ 156 ]), &(nmpcWorkspace.QDy[ 168 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 2366 ]), &(nmpcWorkspace.Dy[ 169 ]), &(nmpcWorkspace.QDy[ 182 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 2548 ]), &(nmpcWorkspace.Dy[ 182 ]), &(nmpcWorkspace.QDy[ 196 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 2730 ]), &(nmpcWorkspace.Dy[ 195 ]), &(nmpcWorkspace.QDy[ 210 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 2912 ]), &(nmpcWorkspace.Dy[ 208 ]), &(nmpcWorkspace.QDy[ 224 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 3094 ]), &(nmpcWorkspace.Dy[ 221 ]), &(nmpcWorkspace.QDy[ 238 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 3276 ]), &(nmpcWorkspace.Dy[ 234 ]), &(nmpcWorkspace.QDy[ 252 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 3458 ]), &(nmpcWorkspace.Dy[ 247 ]), &(nmpcWorkspace.QDy[ 266 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 3640 ]), &(nmpcWorkspace.Dy[ 260 ]), &(nmpcWorkspace.QDy[ 280 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 3822 ]), &(nmpcWorkspace.Dy[ 273 ]), &(nmpcWorkspace.QDy[ 294 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 4004 ]), &(nmpcWorkspace.Dy[ 286 ]), &(nmpcWorkspace.QDy[ 308 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 4186 ]), &(nmpcWorkspace.Dy[ 299 ]), &(nmpcWorkspace.QDy[ 322 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 4368 ]), &(nmpcWorkspace.Dy[ 312 ]), &(nmpcWorkspace.QDy[ 336 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 4550 ]), &(nmpcWorkspace.Dy[ 325 ]), &(nmpcWorkspace.QDy[ 350 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 4732 ]), &(nmpcWorkspace.Dy[ 338 ]), &(nmpcWorkspace.QDy[ 364 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 4914 ]), &(nmpcWorkspace.Dy[ 351 ]), &(nmpcWorkspace.QDy[ 378 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 5096 ]), &(nmpcWorkspace.Dy[ 364 ]), &(nmpcWorkspace.QDy[ 392 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 5278 ]), &(nmpcWorkspace.Dy[ 377 ]), &(nmpcWorkspace.QDy[ 406 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 5460 ]), &(nmpcWorkspace.Dy[ 390 ]), &(nmpcWorkspace.QDy[ 420 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 5642 ]), &(nmpcWorkspace.Dy[ 403 ]), &(nmpcWorkspace.QDy[ 434 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 5824 ]), &(nmpcWorkspace.Dy[ 416 ]), &(nmpcWorkspace.QDy[ 448 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 6006 ]), &(nmpcWorkspace.Dy[ 429 ]), &(nmpcWorkspace.QDy[ 462 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 6188 ]), &(nmpcWorkspace.Dy[ 442 ]), &(nmpcWorkspace.QDy[ 476 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 6370 ]), &(nmpcWorkspace.Dy[ 455 ]), &(nmpcWorkspace.QDy[ 490 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 6552 ]), &(nmpcWorkspace.Dy[ 468 ]), &(nmpcWorkspace.QDy[ 504 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 6734 ]), &(nmpcWorkspace.Dy[ 481 ]), &(nmpcWorkspace.QDy[ 518 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 6916 ]), &(nmpcWorkspace.Dy[ 494 ]), &(nmpcWorkspace.QDy[ 532 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 7098 ]), &(nmpcWorkspace.Dy[ 507 ]), &(nmpcWorkspace.QDy[ 546 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 7280 ]), &(nmpcWorkspace.Dy[ 520 ]), &(nmpcWorkspace.QDy[ 560 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 7462 ]), &(nmpcWorkspace.Dy[ 533 ]), &(nmpcWorkspace.QDy[ 574 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 7644 ]), &(nmpcWorkspace.Dy[ 546 ]), &(nmpcWorkspace.QDy[ 588 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 7826 ]), &(nmpcWorkspace.Dy[ 559 ]), &(nmpcWorkspace.QDy[ 602 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 8008 ]), &(nmpcWorkspace.Dy[ 572 ]), &(nmpcWorkspace.QDy[ 616 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 8190 ]), &(nmpcWorkspace.Dy[ 585 ]), &(nmpcWorkspace.QDy[ 630 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 8372 ]), &(nmpcWorkspace.Dy[ 598 ]), &(nmpcWorkspace.QDy[ 644 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 8554 ]), &(nmpcWorkspace.Dy[ 611 ]), &(nmpcWorkspace.QDy[ 658 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 8736 ]), &(nmpcWorkspace.Dy[ 624 ]), &(nmpcWorkspace.QDy[ 672 ]) );
nmpc_multQDy( &(nmpcWorkspace.Q2[ 8918 ]), &(nmpcWorkspace.Dy[ 637 ]), &(nmpcWorkspace.QDy[ 686 ]) );

nmpcWorkspace.QDy[700] = + nmpcWorkspace.QN2[0]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[1]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[2]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[3]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[4]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[5]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[6]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[7]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[701] = + nmpcWorkspace.QN2[8]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[9]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[10]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[11]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[12]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[13]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[14]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[15]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[702] = + nmpcWorkspace.QN2[16]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[17]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[18]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[19]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[20]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[21]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[22]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[23]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[703] = + nmpcWorkspace.QN2[24]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[25]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[26]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[27]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[28]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[29]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[30]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[31]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[704] = + nmpcWorkspace.QN2[32]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[33]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[34]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[35]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[36]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[37]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[38]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[39]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[705] = + nmpcWorkspace.QN2[40]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[41]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[42]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[43]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[44]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[45]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[46]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[47]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[706] = + nmpcWorkspace.QN2[48]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[49]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[50]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[51]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[52]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[53]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[54]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[55]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[707] = + nmpcWorkspace.QN2[56]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[57]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[58]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[59]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[60]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[61]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[62]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[63]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[708] = + nmpcWorkspace.QN2[64]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[65]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[66]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[67]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[68]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[69]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[70]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[71]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[709] = + nmpcWorkspace.QN2[72]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[73]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[74]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[75]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[76]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[77]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[78]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[79]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[710] = + nmpcWorkspace.QN2[80]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[81]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[82]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[83]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[84]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[85]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[86]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[87]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[711] = + nmpcWorkspace.QN2[88]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[89]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[90]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[91]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[92]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[93]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[94]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[95]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[712] = + nmpcWorkspace.QN2[96]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[97]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[98]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[99]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[100]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[101]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[102]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[103]*nmpcWorkspace.DyN[7];
nmpcWorkspace.QDy[713] = + nmpcWorkspace.QN2[104]*nmpcWorkspace.DyN[0] + nmpcWorkspace.QN2[105]*nmpcWorkspace.DyN[1] + nmpcWorkspace.QN2[106]*nmpcWorkspace.DyN[2] + nmpcWorkspace.QN2[107]*nmpcWorkspace.DyN[3] + nmpcWorkspace.QN2[108]*nmpcWorkspace.DyN[4] + nmpcWorkspace.QN2[109]*nmpcWorkspace.DyN[5] + nmpcWorkspace.QN2[110]*nmpcWorkspace.DyN[6] + nmpcWorkspace.QN2[111]*nmpcWorkspace.DyN[7];

for (lRun2 = 0; lRun2 < 700; ++lRun2)
nmpcWorkspace.QDy[lRun2 + 14] += nmpcWorkspace.Qd[lRun2];


for (lRun2 = 0; lRun2 < 14; ++lRun2)
{
for (lRun4 = 0; lRun4 < 1; ++lRun4)
{
real_t t = 0.0;
for (lRun5 = 0; lRun5 < 700; ++lRun5)
{
t += + nmpcWorkspace.evGx[(lRun5 * 14) + (lRun2)]*nmpcWorkspace.QDy[(lRun5 + 14) + (lRun4)];
}
nmpcWorkspace.g[(lRun2) + (lRun4)] = t;
}
}


for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
for (lRun2 = lRun1; lRun2 < 50; ++lRun2)
{
lRun3 = (((lRun2 + 1) * (lRun2)) / (2)) + (lRun1);
nmpc_multEQDy( &(nmpcWorkspace.E[ lRun3 * 56 ]), &(nmpcWorkspace.QDy[ lRun2 * 14 + 14 ]), &(nmpcWorkspace.g[ lRun1 * 4 + 14 ]) );
}
}

nmpcWorkspace.lb[0] = nmpcWorkspace.Dx0[0];
nmpcWorkspace.lb[1] = nmpcWorkspace.Dx0[1];
nmpcWorkspace.lb[2] = nmpcWorkspace.Dx0[2];
nmpcWorkspace.lb[3] = nmpcWorkspace.Dx0[3];
nmpcWorkspace.lb[4] = nmpcWorkspace.Dx0[4];
nmpcWorkspace.lb[5] = nmpcWorkspace.Dx0[5];
nmpcWorkspace.lb[6] = nmpcWorkspace.Dx0[6];
nmpcWorkspace.lb[7] = nmpcWorkspace.Dx0[7];
nmpcWorkspace.lb[8] = nmpcWorkspace.Dx0[8];
nmpcWorkspace.lb[9] = nmpcWorkspace.Dx0[9];
nmpcWorkspace.lb[10] = nmpcWorkspace.Dx0[10];
nmpcWorkspace.lb[11] = nmpcWorkspace.Dx0[11];
nmpcWorkspace.lb[12] = nmpcWorkspace.Dx0[12];
nmpcWorkspace.lb[13] = nmpcWorkspace.Dx0[13];
nmpcWorkspace.ub[0] = nmpcWorkspace.Dx0[0];
nmpcWorkspace.ub[1] = nmpcWorkspace.Dx0[1];
nmpcWorkspace.ub[2] = nmpcWorkspace.Dx0[2];
nmpcWorkspace.ub[3] = nmpcWorkspace.Dx0[3];
nmpcWorkspace.ub[4] = nmpcWorkspace.Dx0[4];
nmpcWorkspace.ub[5] = nmpcWorkspace.Dx0[5];
nmpcWorkspace.ub[6] = nmpcWorkspace.Dx0[6];
nmpcWorkspace.ub[7] = nmpcWorkspace.Dx0[7];
nmpcWorkspace.ub[8] = nmpcWorkspace.Dx0[8];
nmpcWorkspace.ub[9] = nmpcWorkspace.Dx0[9];
nmpcWorkspace.ub[10] = nmpcWorkspace.Dx0[10];
nmpcWorkspace.ub[11] = nmpcWorkspace.Dx0[11];
nmpcWorkspace.ub[12] = nmpcWorkspace.Dx0[12];
nmpcWorkspace.ub[13] = nmpcWorkspace.Dx0[13];
}

void nmpc_expand(  )
{
int lRun1;
int lRun2;
int lRun3;
nmpcVariables.x[0] += nmpcWorkspace.x[0];
nmpcVariables.x[1] += nmpcWorkspace.x[1];
nmpcVariables.x[2] += nmpcWorkspace.x[2];
nmpcVariables.x[3] += nmpcWorkspace.x[3];
nmpcVariables.x[4] += nmpcWorkspace.x[4];
nmpcVariables.x[5] += nmpcWorkspace.x[5];
nmpcVariables.x[6] += nmpcWorkspace.x[6];
nmpcVariables.x[7] += nmpcWorkspace.x[7];
nmpcVariables.x[8] += nmpcWorkspace.x[8];
nmpcVariables.x[9] += nmpcWorkspace.x[9];
nmpcVariables.x[10] += nmpcWorkspace.x[10];
nmpcVariables.x[11] += nmpcWorkspace.x[11];
nmpcVariables.x[12] += nmpcWorkspace.x[12];
nmpcVariables.x[13] += nmpcWorkspace.x[13];

for (lRun1 = 0; lRun1 < 200; ++lRun1)
nmpcVariables.u[lRun1] += nmpcWorkspace.x[lRun1 + 14];


for (lRun1 = 0; lRun1 < 700; ++lRun1)
{
for (lRun2 = 0; lRun2 < 1; ++lRun2)
{
real_t t = 0.0;
for (lRun3 = 0; lRun3 < 14; ++lRun3)
{
t += + nmpcWorkspace.evGx[(lRun1 * 14) + (lRun3)]*nmpcWorkspace.x[(lRun3) + (lRun2)];
}
nmpcVariables.x[(lRun1 + 14) + (lRun2)] += t + nmpcWorkspace.d[(lRun1) + (lRun2)];
}
}

for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
for (lRun2 = 0; lRun2 < lRun1 + 1; ++lRun2)
{
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
nmpc_multEDu( &(nmpcWorkspace.E[ lRun3 * 56 ]), &(nmpcWorkspace.x[ lRun2 * 4 + 14 ]), &(nmpcVariables.x[ lRun1 * 14 + 14 ]) );
}
}
}

int nmpc_preparationStep(  )
{
int ret;

ret = nmpc_modelSimulation();
nmpc_evaluateObjective(  );
nmpc_condensePrep(  );
return ret;
}

int nmpc_feedbackStep(  )
{
int tmp;

nmpc_condenseFdb(  );

tmp = nmpc_solve( );

nmpc_expand(  );
return tmp;
}

int nmpc_initializeSolver(  )
{
int ret;

/* This is a function which must be called once before any other function call! */


ret = 0;

memset(&nmpcWorkspace, 0, sizeof( nmpcWorkspace ));
return ret;
}

void nmpc_initializeNodesByForwardSimulation(  )
{
int index;
for (index = 0; index < 50; ++index)
{
nmpcWorkspace.state[0] = nmpcVariables.x[index * 14];
nmpcWorkspace.state[1] = nmpcVariables.x[index * 14 + 1];
nmpcWorkspace.state[2] = nmpcVariables.x[index * 14 + 2];
nmpcWorkspace.state[3] = nmpcVariables.x[index * 14 + 3];
nmpcWorkspace.state[4] = nmpcVariables.x[index * 14 + 4];
nmpcWorkspace.state[5] = nmpcVariables.x[index * 14 + 5];
nmpcWorkspace.state[6] = nmpcVariables.x[index * 14 + 6];
nmpcWorkspace.state[7] = nmpcVariables.x[index * 14 + 7];
nmpcWorkspace.state[8] = nmpcVariables.x[index * 14 + 8];
nmpcWorkspace.state[9] = nmpcVariables.x[index * 14 + 9];
nmpcWorkspace.state[10] = nmpcVariables.x[index * 14 + 10];
nmpcWorkspace.state[11] = nmpcVariables.x[index * 14 + 11];
nmpcWorkspace.state[12] = nmpcVariables.x[index * 14 + 12];
nmpcWorkspace.state[13] = nmpcVariables.x[index * 14 + 13];
nmpcWorkspace.state[266] = nmpcVariables.u[index * 4];
nmpcWorkspace.state[267] = nmpcVariables.u[index * 4 + 1];
nmpcWorkspace.state[268] = nmpcVariables.u[index * 4 + 2];
nmpcWorkspace.state[269] = nmpcVariables.u[index * 4 + 3];
nmpcWorkspace.state[270] = nmpcVariables.od[index * 6];
nmpcWorkspace.state[271] = nmpcVariables.od[index * 6 + 1];
nmpcWorkspace.state[272] = nmpcVariables.od[index * 6 + 2];
nmpcWorkspace.state[273] = nmpcVariables.od[index * 6 + 3];
nmpcWorkspace.state[274] = nmpcVariables.od[index * 6 + 4];
nmpcWorkspace.state[275] = nmpcVariables.od[index * 6 + 5];

nmpc_integrate(nmpcWorkspace.state, index == 0);

nmpcVariables.x[index * 14 + 14] = nmpcWorkspace.state[0];
nmpcVariables.x[index * 14 + 15] = nmpcWorkspace.state[1];
nmpcVariables.x[index * 14 + 16] = nmpcWorkspace.state[2];
nmpcVariables.x[index * 14 + 17] = nmpcWorkspace.state[3];
nmpcVariables.x[index * 14 + 18] = nmpcWorkspace.state[4];
nmpcVariables.x[index * 14 + 19] = nmpcWorkspace.state[5];
nmpcVariables.x[index * 14 + 20] = nmpcWorkspace.state[6];
nmpcVariables.x[index * 14 + 21] = nmpcWorkspace.state[7];
nmpcVariables.x[index * 14 + 22] = nmpcWorkspace.state[8];
nmpcVariables.x[index * 14 + 23] = nmpcWorkspace.state[9];
nmpcVariables.x[index * 14 + 24] = nmpcWorkspace.state[10];
nmpcVariables.x[index * 14 + 25] = nmpcWorkspace.state[11];
nmpcVariables.x[index * 14 + 26] = nmpcWorkspace.state[12];
nmpcVariables.x[index * 14 + 27] = nmpcWorkspace.state[13];
}
}

void nmpc_shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd )
{
int index;
for (index = 0; index < 50; ++index)
{
nmpcVariables.x[index * 14] = nmpcVariables.x[index * 14 + 14];
nmpcVariables.x[index * 14 + 1] = nmpcVariables.x[index * 14 + 15];
nmpcVariables.x[index * 14 + 2] = nmpcVariables.x[index * 14 + 16];
nmpcVariables.x[index * 14 + 3] = nmpcVariables.x[index * 14 + 17];
nmpcVariables.x[index * 14 + 4] = nmpcVariables.x[index * 14 + 18];
nmpcVariables.x[index * 14 + 5] = nmpcVariables.x[index * 14 + 19];
nmpcVariables.x[index * 14 + 6] = nmpcVariables.x[index * 14 + 20];
nmpcVariables.x[index * 14 + 7] = nmpcVariables.x[index * 14 + 21];
nmpcVariables.x[index * 14 + 8] = nmpcVariables.x[index * 14 + 22];
nmpcVariables.x[index * 14 + 9] = nmpcVariables.x[index * 14 + 23];
nmpcVariables.x[index * 14 + 10] = nmpcVariables.x[index * 14 + 24];
nmpcVariables.x[index * 14 + 11] = nmpcVariables.x[index * 14 + 25];
nmpcVariables.x[index * 14 + 12] = nmpcVariables.x[index * 14 + 26];
nmpcVariables.x[index * 14 + 13] = nmpcVariables.x[index * 14 + 27];
}

if (strategy == 1 && xEnd != 0)
{
nmpcVariables.x[700] = xEnd[0];
nmpcVariables.x[701] = xEnd[1];
nmpcVariables.x[702] = xEnd[2];
nmpcVariables.x[703] = xEnd[3];
nmpcVariables.x[704] = xEnd[4];
nmpcVariables.x[705] = xEnd[5];
nmpcVariables.x[706] = xEnd[6];
nmpcVariables.x[707] = xEnd[7];
nmpcVariables.x[708] = xEnd[8];
nmpcVariables.x[709] = xEnd[9];
nmpcVariables.x[710] = xEnd[10];
nmpcVariables.x[711] = xEnd[11];
nmpcVariables.x[712] = xEnd[12];
nmpcVariables.x[713] = xEnd[13];
}
else if (strategy == 2) 
{
nmpcWorkspace.state[0] = nmpcVariables.x[700];
nmpcWorkspace.state[1] = nmpcVariables.x[701];
nmpcWorkspace.state[2] = nmpcVariables.x[702];
nmpcWorkspace.state[3] = nmpcVariables.x[703];
nmpcWorkspace.state[4] = nmpcVariables.x[704];
nmpcWorkspace.state[5] = nmpcVariables.x[705];
nmpcWorkspace.state[6] = nmpcVariables.x[706];
nmpcWorkspace.state[7] = nmpcVariables.x[707];
nmpcWorkspace.state[8] = nmpcVariables.x[708];
nmpcWorkspace.state[9] = nmpcVariables.x[709];
nmpcWorkspace.state[10] = nmpcVariables.x[710];
nmpcWorkspace.state[11] = nmpcVariables.x[711];
nmpcWorkspace.state[12] = nmpcVariables.x[712];
nmpcWorkspace.state[13] = nmpcVariables.x[713];
if (uEnd != 0)
{
nmpcWorkspace.state[266] = uEnd[0];
nmpcWorkspace.state[267] = uEnd[1];
nmpcWorkspace.state[268] = uEnd[2];
nmpcWorkspace.state[269] = uEnd[3];
}
else
{
nmpcWorkspace.state[266] = nmpcVariables.u[196];
nmpcWorkspace.state[267] = nmpcVariables.u[197];
nmpcWorkspace.state[268] = nmpcVariables.u[198];
nmpcWorkspace.state[269] = nmpcVariables.u[199];
}
nmpcWorkspace.state[270] = nmpcVariables.od[300];
nmpcWorkspace.state[271] = nmpcVariables.od[301];
nmpcWorkspace.state[272] = nmpcVariables.od[302];
nmpcWorkspace.state[273] = nmpcVariables.od[303];
nmpcWorkspace.state[274] = nmpcVariables.od[304];
nmpcWorkspace.state[275] = nmpcVariables.od[305];

nmpc_integrate(nmpcWorkspace.state, 1);

nmpcVariables.x[700] = nmpcWorkspace.state[0];
nmpcVariables.x[701] = nmpcWorkspace.state[1];
nmpcVariables.x[702] = nmpcWorkspace.state[2];
nmpcVariables.x[703] = nmpcWorkspace.state[3];
nmpcVariables.x[704] = nmpcWorkspace.state[4];
nmpcVariables.x[705] = nmpcWorkspace.state[5];
nmpcVariables.x[706] = nmpcWorkspace.state[6];
nmpcVariables.x[707] = nmpcWorkspace.state[7];
nmpcVariables.x[708] = nmpcWorkspace.state[8];
nmpcVariables.x[709] = nmpcWorkspace.state[9];
nmpcVariables.x[710] = nmpcWorkspace.state[10];
nmpcVariables.x[711] = nmpcWorkspace.state[11];
nmpcVariables.x[712] = nmpcWorkspace.state[12];
nmpcVariables.x[713] = nmpcWorkspace.state[13];
}
}

void nmpc_shiftControls( real_t* const uEnd )
{
int index;
for (index = 0; index < 49; ++index)
{
nmpcVariables.u[index * 4] = nmpcVariables.u[index * 4 + 4];
nmpcVariables.u[index * 4 + 1] = nmpcVariables.u[index * 4 + 5];
nmpcVariables.u[index * 4 + 2] = nmpcVariables.u[index * 4 + 6];
nmpcVariables.u[index * 4 + 3] = nmpcVariables.u[index * 4 + 7];
}

if (uEnd != 0)
{
nmpcVariables.u[196] = uEnd[0];
nmpcVariables.u[197] = uEnd[1];
nmpcVariables.u[198] = uEnd[2];
nmpcVariables.u[199] = uEnd[3];
}
}

real_t nmpc_getKKT(  )
{
real_t kkt;

int index;
real_t prd;

kkt = + nmpcWorkspace.g[0]*nmpcWorkspace.x[0] + nmpcWorkspace.g[1]*nmpcWorkspace.x[1] + nmpcWorkspace.g[2]*nmpcWorkspace.x[2] + nmpcWorkspace.g[3]*nmpcWorkspace.x[3] + nmpcWorkspace.g[4]*nmpcWorkspace.x[4] + nmpcWorkspace.g[5]*nmpcWorkspace.x[5] + nmpcWorkspace.g[6]*nmpcWorkspace.x[6] + nmpcWorkspace.g[7]*nmpcWorkspace.x[7] + nmpcWorkspace.g[8]*nmpcWorkspace.x[8] + nmpcWorkspace.g[9]*nmpcWorkspace.x[9] + nmpcWorkspace.g[10]*nmpcWorkspace.x[10] + nmpcWorkspace.g[11]*nmpcWorkspace.x[11] + nmpcWorkspace.g[12]*nmpcWorkspace.x[12] + nmpcWorkspace.g[13]*nmpcWorkspace.x[13] + nmpcWorkspace.g[14]*nmpcWorkspace.x[14] + nmpcWorkspace.g[15]*nmpcWorkspace.x[15] + nmpcWorkspace.g[16]*nmpcWorkspace.x[16] + nmpcWorkspace.g[17]*nmpcWorkspace.x[17] + nmpcWorkspace.g[18]*nmpcWorkspace.x[18] + nmpcWorkspace.g[19]*nmpcWorkspace.x[19] + nmpcWorkspace.g[20]*nmpcWorkspace.x[20] + nmpcWorkspace.g[21]*nmpcWorkspace.x[21] + nmpcWorkspace.g[22]*nmpcWorkspace.x[22] + nmpcWorkspace.g[23]*nmpcWorkspace.x[23] + nmpcWorkspace.g[24]*nmpcWorkspace.x[24] + nmpcWorkspace.g[25]*nmpcWorkspace.x[25] + nmpcWorkspace.g[26]*nmpcWorkspace.x[26] + nmpcWorkspace.g[27]*nmpcWorkspace.x[27] + nmpcWorkspace.g[28]*nmpcWorkspace.x[28] + nmpcWorkspace.g[29]*nmpcWorkspace.x[29] + nmpcWorkspace.g[30]*nmpcWorkspace.x[30] + nmpcWorkspace.g[31]*nmpcWorkspace.x[31] + nmpcWorkspace.g[32]*nmpcWorkspace.x[32] + nmpcWorkspace.g[33]*nmpcWorkspace.x[33] + nmpcWorkspace.g[34]*nmpcWorkspace.x[34] + nmpcWorkspace.g[35]*nmpcWorkspace.x[35] + nmpcWorkspace.g[36]*nmpcWorkspace.x[36] + nmpcWorkspace.g[37]*nmpcWorkspace.x[37] + nmpcWorkspace.g[38]*nmpcWorkspace.x[38] + nmpcWorkspace.g[39]*nmpcWorkspace.x[39] + nmpcWorkspace.g[40]*nmpcWorkspace.x[40] + nmpcWorkspace.g[41]*nmpcWorkspace.x[41] + nmpcWorkspace.g[42]*nmpcWorkspace.x[42] + nmpcWorkspace.g[43]*nmpcWorkspace.x[43] + nmpcWorkspace.g[44]*nmpcWorkspace.x[44] + nmpcWorkspace.g[45]*nmpcWorkspace.x[45] + nmpcWorkspace.g[46]*nmpcWorkspace.x[46] + nmpcWorkspace.g[47]*nmpcWorkspace.x[47] + nmpcWorkspace.g[48]*nmpcWorkspace.x[48] + nmpcWorkspace.g[49]*nmpcWorkspace.x[49] + nmpcWorkspace.g[50]*nmpcWorkspace.x[50] + nmpcWorkspace.g[51]*nmpcWorkspace.x[51] + nmpcWorkspace.g[52]*nmpcWorkspace.x[52] + nmpcWorkspace.g[53]*nmpcWorkspace.x[53] + nmpcWorkspace.g[54]*nmpcWorkspace.x[54] + nmpcWorkspace.g[55]*nmpcWorkspace.x[55] + nmpcWorkspace.g[56]*nmpcWorkspace.x[56] + nmpcWorkspace.g[57]*nmpcWorkspace.x[57] + nmpcWorkspace.g[58]*nmpcWorkspace.x[58] + nmpcWorkspace.g[59]*nmpcWorkspace.x[59] + nmpcWorkspace.g[60]*nmpcWorkspace.x[60] + nmpcWorkspace.g[61]*nmpcWorkspace.x[61] + nmpcWorkspace.g[62]*nmpcWorkspace.x[62] + nmpcWorkspace.g[63]*nmpcWorkspace.x[63] + nmpcWorkspace.g[64]*nmpcWorkspace.x[64] + nmpcWorkspace.g[65]*nmpcWorkspace.x[65] + nmpcWorkspace.g[66]*nmpcWorkspace.x[66] + nmpcWorkspace.g[67]*nmpcWorkspace.x[67] + nmpcWorkspace.g[68]*nmpcWorkspace.x[68] + nmpcWorkspace.g[69]*nmpcWorkspace.x[69] + nmpcWorkspace.g[70]*nmpcWorkspace.x[70] + nmpcWorkspace.g[71]*nmpcWorkspace.x[71] + nmpcWorkspace.g[72]*nmpcWorkspace.x[72] + nmpcWorkspace.g[73]*nmpcWorkspace.x[73] + nmpcWorkspace.g[74]*nmpcWorkspace.x[74] + nmpcWorkspace.g[75]*nmpcWorkspace.x[75] + nmpcWorkspace.g[76]*nmpcWorkspace.x[76] + nmpcWorkspace.g[77]*nmpcWorkspace.x[77] + nmpcWorkspace.g[78]*nmpcWorkspace.x[78] + nmpcWorkspace.g[79]*nmpcWorkspace.x[79] + nmpcWorkspace.g[80]*nmpcWorkspace.x[80] + nmpcWorkspace.g[81]*nmpcWorkspace.x[81] + nmpcWorkspace.g[82]*nmpcWorkspace.x[82] + nmpcWorkspace.g[83]*nmpcWorkspace.x[83] + nmpcWorkspace.g[84]*nmpcWorkspace.x[84] + nmpcWorkspace.g[85]*nmpcWorkspace.x[85] + nmpcWorkspace.g[86]*nmpcWorkspace.x[86] + nmpcWorkspace.g[87]*nmpcWorkspace.x[87] + nmpcWorkspace.g[88]*nmpcWorkspace.x[88] + nmpcWorkspace.g[89]*nmpcWorkspace.x[89] + nmpcWorkspace.g[90]*nmpcWorkspace.x[90] + nmpcWorkspace.g[91]*nmpcWorkspace.x[91] + nmpcWorkspace.g[92]*nmpcWorkspace.x[92] + nmpcWorkspace.g[93]*nmpcWorkspace.x[93] + nmpcWorkspace.g[94]*nmpcWorkspace.x[94] + nmpcWorkspace.g[95]*nmpcWorkspace.x[95] + nmpcWorkspace.g[96]*nmpcWorkspace.x[96] + nmpcWorkspace.g[97]*nmpcWorkspace.x[97] + nmpcWorkspace.g[98]*nmpcWorkspace.x[98] + nmpcWorkspace.g[99]*nmpcWorkspace.x[99] + nmpcWorkspace.g[100]*nmpcWorkspace.x[100] + nmpcWorkspace.g[101]*nmpcWorkspace.x[101] + nmpcWorkspace.g[102]*nmpcWorkspace.x[102] + nmpcWorkspace.g[103]*nmpcWorkspace.x[103] + nmpcWorkspace.g[104]*nmpcWorkspace.x[104] + nmpcWorkspace.g[105]*nmpcWorkspace.x[105] + nmpcWorkspace.g[106]*nmpcWorkspace.x[106] + nmpcWorkspace.g[107]*nmpcWorkspace.x[107] + nmpcWorkspace.g[108]*nmpcWorkspace.x[108] + nmpcWorkspace.g[109]*nmpcWorkspace.x[109] + nmpcWorkspace.g[110]*nmpcWorkspace.x[110] + nmpcWorkspace.g[111]*nmpcWorkspace.x[111] + nmpcWorkspace.g[112]*nmpcWorkspace.x[112] + nmpcWorkspace.g[113]*nmpcWorkspace.x[113] + nmpcWorkspace.g[114]*nmpcWorkspace.x[114] + nmpcWorkspace.g[115]*nmpcWorkspace.x[115] + nmpcWorkspace.g[116]*nmpcWorkspace.x[116] + nmpcWorkspace.g[117]*nmpcWorkspace.x[117] + nmpcWorkspace.g[118]*nmpcWorkspace.x[118] + nmpcWorkspace.g[119]*nmpcWorkspace.x[119] + nmpcWorkspace.g[120]*nmpcWorkspace.x[120] + nmpcWorkspace.g[121]*nmpcWorkspace.x[121] + nmpcWorkspace.g[122]*nmpcWorkspace.x[122] + nmpcWorkspace.g[123]*nmpcWorkspace.x[123] + nmpcWorkspace.g[124]*nmpcWorkspace.x[124] + nmpcWorkspace.g[125]*nmpcWorkspace.x[125] + nmpcWorkspace.g[126]*nmpcWorkspace.x[126] + nmpcWorkspace.g[127]*nmpcWorkspace.x[127] + nmpcWorkspace.g[128]*nmpcWorkspace.x[128] + nmpcWorkspace.g[129]*nmpcWorkspace.x[129] + nmpcWorkspace.g[130]*nmpcWorkspace.x[130] + nmpcWorkspace.g[131]*nmpcWorkspace.x[131] + nmpcWorkspace.g[132]*nmpcWorkspace.x[132] + nmpcWorkspace.g[133]*nmpcWorkspace.x[133] + nmpcWorkspace.g[134]*nmpcWorkspace.x[134] + nmpcWorkspace.g[135]*nmpcWorkspace.x[135] + nmpcWorkspace.g[136]*nmpcWorkspace.x[136] + nmpcWorkspace.g[137]*nmpcWorkspace.x[137] + nmpcWorkspace.g[138]*nmpcWorkspace.x[138] + nmpcWorkspace.g[139]*nmpcWorkspace.x[139] + nmpcWorkspace.g[140]*nmpcWorkspace.x[140] + nmpcWorkspace.g[141]*nmpcWorkspace.x[141] + nmpcWorkspace.g[142]*nmpcWorkspace.x[142] + nmpcWorkspace.g[143]*nmpcWorkspace.x[143] + nmpcWorkspace.g[144]*nmpcWorkspace.x[144] + nmpcWorkspace.g[145]*nmpcWorkspace.x[145] + nmpcWorkspace.g[146]*nmpcWorkspace.x[146] + nmpcWorkspace.g[147]*nmpcWorkspace.x[147] + nmpcWorkspace.g[148]*nmpcWorkspace.x[148] + nmpcWorkspace.g[149]*nmpcWorkspace.x[149] + nmpcWorkspace.g[150]*nmpcWorkspace.x[150] + nmpcWorkspace.g[151]*nmpcWorkspace.x[151] + nmpcWorkspace.g[152]*nmpcWorkspace.x[152] + nmpcWorkspace.g[153]*nmpcWorkspace.x[153] + nmpcWorkspace.g[154]*nmpcWorkspace.x[154] + nmpcWorkspace.g[155]*nmpcWorkspace.x[155] + nmpcWorkspace.g[156]*nmpcWorkspace.x[156] + nmpcWorkspace.g[157]*nmpcWorkspace.x[157] + nmpcWorkspace.g[158]*nmpcWorkspace.x[158] + nmpcWorkspace.g[159]*nmpcWorkspace.x[159] + nmpcWorkspace.g[160]*nmpcWorkspace.x[160] + nmpcWorkspace.g[161]*nmpcWorkspace.x[161] + nmpcWorkspace.g[162]*nmpcWorkspace.x[162] + nmpcWorkspace.g[163]*nmpcWorkspace.x[163] + nmpcWorkspace.g[164]*nmpcWorkspace.x[164] + nmpcWorkspace.g[165]*nmpcWorkspace.x[165] + nmpcWorkspace.g[166]*nmpcWorkspace.x[166] + nmpcWorkspace.g[167]*nmpcWorkspace.x[167] + nmpcWorkspace.g[168]*nmpcWorkspace.x[168] + nmpcWorkspace.g[169]*nmpcWorkspace.x[169] + nmpcWorkspace.g[170]*nmpcWorkspace.x[170] + nmpcWorkspace.g[171]*nmpcWorkspace.x[171] + nmpcWorkspace.g[172]*nmpcWorkspace.x[172] + nmpcWorkspace.g[173]*nmpcWorkspace.x[173] + nmpcWorkspace.g[174]*nmpcWorkspace.x[174] + nmpcWorkspace.g[175]*nmpcWorkspace.x[175] + nmpcWorkspace.g[176]*nmpcWorkspace.x[176] + nmpcWorkspace.g[177]*nmpcWorkspace.x[177] + nmpcWorkspace.g[178]*nmpcWorkspace.x[178] + nmpcWorkspace.g[179]*nmpcWorkspace.x[179] + nmpcWorkspace.g[180]*nmpcWorkspace.x[180] + nmpcWorkspace.g[181]*nmpcWorkspace.x[181] + nmpcWorkspace.g[182]*nmpcWorkspace.x[182] + nmpcWorkspace.g[183]*nmpcWorkspace.x[183] + nmpcWorkspace.g[184]*nmpcWorkspace.x[184] + nmpcWorkspace.g[185]*nmpcWorkspace.x[185] + nmpcWorkspace.g[186]*nmpcWorkspace.x[186] + nmpcWorkspace.g[187]*nmpcWorkspace.x[187] + nmpcWorkspace.g[188]*nmpcWorkspace.x[188] + nmpcWorkspace.g[189]*nmpcWorkspace.x[189] + nmpcWorkspace.g[190]*nmpcWorkspace.x[190] + nmpcWorkspace.g[191]*nmpcWorkspace.x[191] + nmpcWorkspace.g[192]*nmpcWorkspace.x[192] + nmpcWorkspace.g[193]*nmpcWorkspace.x[193] + nmpcWorkspace.g[194]*nmpcWorkspace.x[194] + nmpcWorkspace.g[195]*nmpcWorkspace.x[195] + nmpcWorkspace.g[196]*nmpcWorkspace.x[196] + nmpcWorkspace.g[197]*nmpcWorkspace.x[197] + nmpcWorkspace.g[198]*nmpcWorkspace.x[198] + nmpcWorkspace.g[199]*nmpcWorkspace.x[199] + nmpcWorkspace.g[200]*nmpcWorkspace.x[200] + nmpcWorkspace.g[201]*nmpcWorkspace.x[201] + nmpcWorkspace.g[202]*nmpcWorkspace.x[202] + nmpcWorkspace.g[203]*nmpcWorkspace.x[203] + nmpcWorkspace.g[204]*nmpcWorkspace.x[204] + nmpcWorkspace.g[205]*nmpcWorkspace.x[205] + nmpcWorkspace.g[206]*nmpcWorkspace.x[206] + nmpcWorkspace.g[207]*nmpcWorkspace.x[207] + nmpcWorkspace.g[208]*nmpcWorkspace.x[208] + nmpcWorkspace.g[209]*nmpcWorkspace.x[209] + nmpcWorkspace.g[210]*nmpcWorkspace.x[210] + nmpcWorkspace.g[211]*nmpcWorkspace.x[211] + nmpcWorkspace.g[212]*nmpcWorkspace.x[212] + nmpcWorkspace.g[213]*nmpcWorkspace.x[213];
kkt = fabs( kkt );
for (index = 0; index < 214; ++index)
{
prd = nmpcWorkspace.y[index];
if (prd > 1e-12)
kkt += fabs(nmpcWorkspace.lb[index] * prd);
else if (prd < -1e-12)
kkt += fabs(nmpcWorkspace.ub[index] * prd);
}
return kkt;
}

real_t nmpc_getObjective(  )
{
real_t objVal;

int lRun1;
/** Row vector of size: 13 */
real_t tmpDy[ 13 ];

/** Row vector of size: 8 */
real_t tmpDyN[ 8 ];

for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
nmpcWorkspace.objValueIn[0] = nmpcVariables.x[lRun1 * 14];
nmpcWorkspace.objValueIn[1] = nmpcVariables.x[lRun1 * 14 + 1];
nmpcWorkspace.objValueIn[2] = nmpcVariables.x[lRun1 * 14 + 2];
nmpcWorkspace.objValueIn[3] = nmpcVariables.x[lRun1 * 14 + 3];
nmpcWorkspace.objValueIn[4] = nmpcVariables.x[lRun1 * 14 + 4];
nmpcWorkspace.objValueIn[5] = nmpcVariables.x[lRun1 * 14 + 5];
nmpcWorkspace.objValueIn[6] = nmpcVariables.x[lRun1 * 14 + 6];
nmpcWorkspace.objValueIn[7] = nmpcVariables.x[lRun1 * 14 + 7];
nmpcWorkspace.objValueIn[8] = nmpcVariables.x[lRun1 * 14 + 8];
nmpcWorkspace.objValueIn[9] = nmpcVariables.x[lRun1 * 14 + 9];
nmpcWorkspace.objValueIn[10] = nmpcVariables.x[lRun1 * 14 + 10];
nmpcWorkspace.objValueIn[11] = nmpcVariables.x[lRun1 * 14 + 11];
nmpcWorkspace.objValueIn[12] = nmpcVariables.x[lRun1 * 14 + 12];
nmpcWorkspace.objValueIn[13] = nmpcVariables.x[lRun1 * 14 + 13];
nmpcWorkspace.objValueIn[14] = nmpcVariables.u[lRun1 * 4];
nmpcWorkspace.objValueIn[15] = nmpcVariables.u[lRun1 * 4 + 1];
nmpcWorkspace.objValueIn[16] = nmpcVariables.u[lRun1 * 4 + 2];
nmpcWorkspace.objValueIn[17] = nmpcVariables.u[lRun1 * 4 + 3];
nmpcWorkspace.objValueIn[18] = nmpcVariables.od[lRun1 * 6];
nmpcWorkspace.objValueIn[19] = nmpcVariables.od[lRun1 * 6 + 1];
nmpcWorkspace.objValueIn[20] = nmpcVariables.od[lRun1 * 6 + 2];
nmpcWorkspace.objValueIn[21] = nmpcVariables.od[lRun1 * 6 + 3];
nmpcWorkspace.objValueIn[22] = nmpcVariables.od[lRun1 * 6 + 4];
nmpcWorkspace.objValueIn[23] = nmpcVariables.od[lRun1 * 6 + 5];

nmpc_evaluateLSQ( nmpcWorkspace.objValueIn, nmpcWorkspace.objValueOut );
nmpcWorkspace.Dy[lRun1 * 13] = nmpcWorkspace.objValueOut[0] - nmpcVariables.y[lRun1 * 13];
nmpcWorkspace.Dy[lRun1 * 13 + 1] = nmpcWorkspace.objValueOut[1] - nmpcVariables.y[lRun1 * 13 + 1];
nmpcWorkspace.Dy[lRun1 * 13 + 2] = nmpcWorkspace.objValueOut[2] - nmpcVariables.y[lRun1 * 13 + 2];
nmpcWorkspace.Dy[lRun1 * 13 + 3] = nmpcWorkspace.objValueOut[3] - nmpcVariables.y[lRun1 * 13 + 3];
nmpcWorkspace.Dy[lRun1 * 13 + 4] = nmpcWorkspace.objValueOut[4] - nmpcVariables.y[lRun1 * 13 + 4];
nmpcWorkspace.Dy[lRun1 * 13 + 5] = nmpcWorkspace.objValueOut[5] - nmpcVariables.y[lRun1 * 13 + 5];
nmpcWorkspace.Dy[lRun1 * 13 + 6] = nmpcWorkspace.objValueOut[6] - nmpcVariables.y[lRun1 * 13 + 6];
nmpcWorkspace.Dy[lRun1 * 13 + 7] = nmpcWorkspace.objValueOut[7] - nmpcVariables.y[lRun1 * 13 + 7];
nmpcWorkspace.Dy[lRun1 * 13 + 8] = nmpcWorkspace.objValueOut[8] - nmpcVariables.y[lRun1 * 13 + 8];
nmpcWorkspace.Dy[lRun1 * 13 + 9] = nmpcWorkspace.objValueOut[9] - nmpcVariables.y[lRun1 * 13 + 9];
nmpcWorkspace.Dy[lRun1 * 13 + 10] = nmpcWorkspace.objValueOut[10] - nmpcVariables.y[lRun1 * 13 + 10];
nmpcWorkspace.Dy[lRun1 * 13 + 11] = nmpcWorkspace.objValueOut[11] - nmpcVariables.y[lRun1 * 13 + 11];
nmpcWorkspace.Dy[lRun1 * 13 + 12] = nmpcWorkspace.objValueOut[12] - nmpcVariables.y[lRun1 * 13 + 12];
}
nmpcWorkspace.objValueIn[0] = nmpcVariables.x[700];
nmpcWorkspace.objValueIn[1] = nmpcVariables.x[701];
nmpcWorkspace.objValueIn[2] = nmpcVariables.x[702];
nmpcWorkspace.objValueIn[3] = nmpcVariables.x[703];
nmpcWorkspace.objValueIn[4] = nmpcVariables.x[704];
nmpcWorkspace.objValueIn[5] = nmpcVariables.x[705];
nmpcWorkspace.objValueIn[6] = nmpcVariables.x[706];
nmpcWorkspace.objValueIn[7] = nmpcVariables.x[707];
nmpcWorkspace.objValueIn[8] = nmpcVariables.x[708];
nmpcWorkspace.objValueIn[9] = nmpcVariables.x[709];
nmpcWorkspace.objValueIn[10] = nmpcVariables.x[710];
nmpcWorkspace.objValueIn[11] = nmpcVariables.x[711];
nmpcWorkspace.objValueIn[12] = nmpcVariables.x[712];
nmpcWorkspace.objValueIn[13] = nmpcVariables.x[713];
nmpcWorkspace.objValueIn[14] = nmpcVariables.od[300];
nmpcWorkspace.objValueIn[15] = nmpcVariables.od[301];
nmpcWorkspace.objValueIn[16] = nmpcVariables.od[302];
nmpcWorkspace.objValueIn[17] = nmpcVariables.od[303];
nmpcWorkspace.objValueIn[18] = nmpcVariables.od[304];
nmpcWorkspace.objValueIn[19] = nmpcVariables.od[305];
nmpc_evaluateLSQEndTerm( nmpcWorkspace.objValueIn, nmpcWorkspace.objValueOut );
nmpcWorkspace.DyN[0] = nmpcWorkspace.objValueOut[0] - nmpcVariables.yN[0];
nmpcWorkspace.DyN[1] = nmpcWorkspace.objValueOut[1] - nmpcVariables.yN[1];
nmpcWorkspace.DyN[2] = nmpcWorkspace.objValueOut[2] - nmpcVariables.yN[2];
nmpcWorkspace.DyN[3] = nmpcWorkspace.objValueOut[3] - nmpcVariables.yN[3];
nmpcWorkspace.DyN[4] = nmpcWorkspace.objValueOut[4] - nmpcVariables.yN[4];
nmpcWorkspace.DyN[5] = nmpcWorkspace.objValueOut[5] - nmpcVariables.yN[5];
nmpcWorkspace.DyN[6] = nmpcWorkspace.objValueOut[6] - nmpcVariables.yN[6];
nmpcWorkspace.DyN[7] = nmpcWorkspace.objValueOut[7] - nmpcVariables.yN[7];
objVal = 0.0000000000000000e+00;
for (lRun1 = 0; lRun1 < 50; ++lRun1)
{
tmpDy[0] = + nmpcWorkspace.Dy[lRun1 * 13]*nmpcVariables.W[0];
tmpDy[1] = + nmpcWorkspace.Dy[lRun1 * 13 + 1]*nmpcVariables.W[14];
tmpDy[2] = + nmpcWorkspace.Dy[lRun1 * 13 + 2]*nmpcVariables.W[28];
tmpDy[3] = + nmpcWorkspace.Dy[lRun1 * 13 + 3]*nmpcVariables.W[42];
tmpDy[4] = + nmpcWorkspace.Dy[lRun1 * 13 + 4]*nmpcVariables.W[56];
tmpDy[5] = + nmpcWorkspace.Dy[lRun1 * 13 + 5]*nmpcVariables.W[70];
tmpDy[6] = + nmpcWorkspace.Dy[lRun1 * 13 + 6]*nmpcVariables.W[84];
tmpDy[7] = + nmpcWorkspace.Dy[lRun1 * 13 + 7]*nmpcVariables.W[98];
tmpDy[8] = + nmpcWorkspace.Dy[lRun1 * 13 + 8]*nmpcVariables.W[112];
tmpDy[9] = + nmpcWorkspace.Dy[lRun1 * 13 + 9]*nmpcVariables.W[126];
tmpDy[10] = + nmpcWorkspace.Dy[lRun1 * 13 + 10]*nmpcVariables.W[140];
tmpDy[11] = + nmpcWorkspace.Dy[lRun1 * 13 + 11]*nmpcVariables.W[154];
tmpDy[12] = + nmpcWorkspace.Dy[lRun1 * 13 + 12]*nmpcVariables.W[168];
objVal += + nmpcWorkspace.Dy[lRun1 * 13]*tmpDy[0] + nmpcWorkspace.Dy[lRun1 * 13 + 1]*tmpDy[1] + nmpcWorkspace.Dy[lRun1 * 13 + 2]*tmpDy[2] + nmpcWorkspace.Dy[lRun1 * 13 + 3]*tmpDy[3] + nmpcWorkspace.Dy[lRun1 * 13 + 4]*tmpDy[4] + nmpcWorkspace.Dy[lRun1 * 13 + 5]*tmpDy[5] + nmpcWorkspace.Dy[lRun1 * 13 + 6]*tmpDy[6] + nmpcWorkspace.Dy[lRun1 * 13 + 7]*tmpDy[7] + nmpcWorkspace.Dy[lRun1 * 13 + 8]*tmpDy[8] + nmpcWorkspace.Dy[lRun1 * 13 + 9]*tmpDy[9] + nmpcWorkspace.Dy[lRun1 * 13 + 10]*tmpDy[10] + nmpcWorkspace.Dy[lRun1 * 13 + 11]*tmpDy[11] + nmpcWorkspace.Dy[lRun1 * 13 + 12]*tmpDy[12];
}

tmpDyN[0] = + nmpcWorkspace.DyN[0]*nmpcVariables.WN[0];
tmpDyN[1] = + nmpcWorkspace.DyN[1]*nmpcVariables.WN[9];
tmpDyN[2] = + nmpcWorkspace.DyN[2]*nmpcVariables.WN[18];
tmpDyN[3] = + nmpcWorkspace.DyN[3]*nmpcVariables.WN[27];
tmpDyN[4] = + nmpcWorkspace.DyN[4]*nmpcVariables.WN[36];
tmpDyN[5] = + nmpcWorkspace.DyN[5]*nmpcVariables.WN[45];
tmpDyN[6] = + nmpcWorkspace.DyN[6]*nmpcVariables.WN[54];
tmpDyN[7] = + nmpcWorkspace.DyN[7]*nmpcVariables.WN[63];
objVal += + nmpcWorkspace.DyN[0]*tmpDyN[0] + nmpcWorkspace.DyN[1]*tmpDyN[1] + nmpcWorkspace.DyN[2]*tmpDyN[2] + nmpcWorkspace.DyN[3]*tmpDyN[3] + nmpcWorkspace.DyN[4]*tmpDyN[4] + nmpcWorkspace.DyN[5]*tmpDyN[5] + nmpcWorkspace.DyN[6]*tmpDyN[6] + nmpcWorkspace.DyN[7]*tmpDyN[7];

objVal *= 0.5;
return objVal;
}

