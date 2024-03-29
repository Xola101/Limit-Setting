# $Id: limit_es.py 513259 2012-08-10 15:41:28Z adye $
# To be included by limit_submit to specify the parameters for
# combined workspaces

#wsfiles     = 'H_comb_%d_AaronArmbruster.root'
#wsfiles     = 'combined_lp2011v2/%d.root'
#wsfiles     = 'combined_lhc3/%d/lhc_atlas.root'
#wsfiles     = 'lp2011v7/lp2_cb/%d.root'
#wsfiles     = 'ccl2011v4/ccl_cb/%g.root'
#wsfiles     = 'paper2011v6/ccl_cb/%g.root'
#wsfiles     = 'moriond2012v4es/mor_cb/%g.root'
#wsfiles     = 'spring2012v2r1es/spr_cb/%g.root'
#wsfiles     = 'spring4r1ic44es/spr_cb/%g.root'
#wsfiles     = 'ICHEP2012/v4/ichep4es/%g.root'
#wsfiles     = 'Discovery2012_AWS/v2_es/%g_combined_lvlv_es.root'
wsfiles     = 'Discovery2012/v3/ES/%g.root'
if "-V" not in opt: dsver += "es"
workspaces="comb"   # add "es" flag to dsver, not here

parm = {
# M      mumin mumax  toys/hour  points
# ===    ===== =====  =========  ======
  110:   [ 0,     0,       8.8, "10/0.2"             ],
  111:   [ 0,     0,       8.8, "10/0.2"             ],
  112:   [ 0,     0,       8.8, "10/0.2"             ],
  113:   [ 0,     0,       8.8, "10/0.2"             ],
  114:   [ 0,     0,       8.8, "10/0.2"             ],
  115:   [ 0,     0,       4.6, "10/0.2"             ],
  116:   [ 0,     0,       4.6, "10/0.2"             ],
  117:   [ 0,     0,       4.6, "10/0.2"             ],
  118:   [ 0,     0,       4.6, "10/0.2"             ],
  119:   [ 0,     0,       4.6, "9/0.2,10"           ],
  120:   [ 0,     0,       8.3, "0.6/0.2,2/0.1,5.6/0.2,6,10/2,/5"           ], # 0.7-5.5
  120.5: [ 0,     0,       8.3, "0.6/0.2,2/0.1,5.6/0.2,6,10/2,/5"           ], # 0.7-5.5
  121:   [ 0,     0,       8.3, "9/0.2,10"           ],
  121.5: [ 0,     0,       8.3, "9/0.2,10"           ],
  122:   [ 0,     0,       8.3, "8/0.2,10/2"         ],
  122.5: [ 0,     0,       8.3, "8/0.2,10/2"         ],
  123:   [ 0,     0,       8.3, "8/0.2,10/2"         ],
  123.5: [ 0,     0,       8.3, "8/0.2,10/2"         ],
  124:   [ 0,     0,       8.3, "7/0.2,6,10/2"       ],
  124.5: [ 0,     0,       8.3, "7/0.2,6,10/2"       ],
  125:   [ 0,     0,       8.3, "2/0.1,5/0.2,6,10/2" ],
  125.5: [ 0,     0,       8.3, "2/0.1,5/0.2,6,10/2" ],
  126:   [ 0,     0,       8.3, "2/0.1,5/0.2,6,10/2" ],
  126.5: [ 0,     0,       8.2, "2/0.1,5/0.2,6,10/2" ],
  127:   [ 0,     0,       8.3, "2/0.1,6/0.2,10/2"   ],
  127.5: [ 0,     0,       8.3, "2/0.1,6/0.2,10/2"   ],
  128:   [ 0,     0,       8.3, "2/0.1,6/0.2,10/2"   ],
  128.5: [ 0,     0,       8.3, "2/0.1,6/0.2,10/2"   ],
  129:   [ 0,     0,       8.3, "2/0.1,5/0.2,6,10/2" ],
  129.5: [ 0,     0,       8.3, "2/0.1,5/0.2,6,10/2" ],
  130:   [ 0,     0,       8.3, "2/0.1,5/0.2,6,10/2" ],
  131:   [ 0,     0,       8.3, "2/0.1,5/0.2,6,10/2" ],
  132:   [ 0,     0,       4.4, "2/0.1,5/0.2,6,10/2" ],
  133:   [ 0,     0,       4.4, "2/0.1,5/0.2,6,10/2" ],
  134:   [ 0,     0,       4.4, "2/0.1,4/0.2,10/2"   ],
  135:   [ 0,     0,       7.3, "2/0.1,5/0.2,6,10/2" ],
  136:   [ 0,     0,       7.3, "2/0.1,4/0.2,10/2"   ],
  137:   [ 0,     0,       7.3, "2/0.1,4/0.2,10/2"   ],
  138:   [ 0,     0,       7.3, "2/0.1,4/0.2,10/2"   ],
  139:   [ 0,     0,       7.3, "2/0.1,4/0.2,10/2"   ],
  140:   [ 0,     0,      13.0, "1.6/0.1,2,6/2"   ], # 0.2-1.6
  141:   [ 0,     0,      13.0, "2/0.1,4/0.2,10/2"   ],
  142:   [ 0,     0,      13.0, "2/0.1,4/0.2,10/2"   ],
  143:   [ 0,     0,      13.0, "2/0.1,4/0.2,10/2"   ],
  144:   [ 0,     0,      13.0, "2/0.1,4/0.2,10/2"   ],
  145:   [ 0,     0,      19.0, "2/0.1,4/0.2,10/2"   ],
  146:   [ 0,     0,      19.0, "2/0.1,4/0.2,10/2"   ],
  147:   [ 0,     0,      19.0, "2/0.1,4/0.2,10/2"   ],
  148:   [ 0,     0,      19.0, "2/0.1,4/0.2,10/2"   ],
  149:   [ 0,     0,      19.0, "2/0.1,4/0.2,10/2"   ],
  150:   [ 0,     0,       6.0, "2/0.1,4/0.2,10/2"   ],
  152:   [ 0,     0,       6.0, "2/0.1,4/0.2,10/2"   ],
  154:   [ 0,     0,       6.0, "2/0.1,4/0.2,10/2"   ],
  156:   [ 0,     0,      53.1, "2/0.1,4/0.2,10/2"   ],
  158:   [ 0,     0,      53.1, "2/0.1,4/0.2,10/2"   ],
  160:   [ 0,     0,      40.5, "0.1,0.2/0.05,1/0.1,2/0.5,3"   ], # 0.1-0.9
  162:   [ 0,     0,      40.5, "2/0.1,4/0.2,10/2"   ],
  164:   [ 0,     0,      40.5, "2/0.1,4/0.2,10/2"   ],
  166:   [ 0,     0,      40.5, "2/0.1,4/0.2,10/2"   ],
  168:   [ 0,     0,      40.5, "2/0.1,4/0.2,10/2"   ],
  170:   [ 0,     0,      50.3, "2/0.1,4/0.2,10/2"   ],
  172:   [ 0,     0,      50.3, "2/0.1,4/0.2,10/2"   ],
  174:   [ 0,     0,      50.3, "2/0.1,10/0.2"       ],
  176:   [ 0,     0,      50.3, "2/0.1,10/0.2"       ],
  178:   [ 0,     0,      50.3, "2/0.1,10/0.2"       ],
  180:   [ 0,     0,      40.3, "2/0.1,4/0.2,10/1"   ],
  182:   [ 0,     0,      40.3, "2/0.1,10/0.2"       ],
  184:   [ 0,     0,      40.3, "2/0.1,10/0.2"       ],
  186:   [ 0,     0,      40.3, "2/0.1,10/0.2"       ],
  188:   [ 0,     0,      40.3, "2/0.1,10/0.2"       ],
  190:   [ 0,     0,      53.7, "2/0.1,10/0.2"       ],
  192:   [ 0,     0,      53.7, "2/0.1,10/0.2"       ],
  194:   [ 0,     0,      53.7, "2/0.1,10/0.2"       ],
  196:   [ 0,     0,      53.7, "2/0.1,10/0.2"       ],
  198:   [ 0,     0,      53.7, "2/0.1,10/0.2"       ],
  200:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  202:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  204:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  206:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  208:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  210:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  212:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  214:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  216:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  218:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  220:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  222:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  224:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  226:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  228:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  230:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  232:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  234:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  236:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  238:   [ 0,     0,       247, "2/0.2,6/0.5,8/2"    ],
  240:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  242:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  244:   [ 0,     0,        77, "0.2,0.4,0.5,1.4/0.1,2/0.2,6/2"  ],
  246:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  248:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  250:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  252:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  254:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  256:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  258:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  260:   [ 0,     0,        77, "0.2,2/0.1,2.4/0.2,4,8/2"    ], # 0.3-2.4
  262:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  264:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  266:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  268:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  270:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  272:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  274:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  276:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  278:   [ 0,     0,        77, "2/0.2,6/0.5,8/2"    ],
  280:   [ 0,     0,        77, "2/0.2,4/0.5,8/1"    ],
  282:   [ 0,     0,        77, "2/0.2,4/0.5,8/1"    ],
  284:   [ 0,     0,        77, "2/0.2,4/0.5,8/1"    ],
  286:   [ 0,     0,        77, "2/0.2,4/0.5,8/1"    ],
  288:   [ 0,     0,        77, "2/0.2,4/0.5,8/1"    ],
  290:   [ 0,     0,        77, "2/0.2,4/0.5,8/1"    ],
  295:   [ 0,     0,        77, "2/0.2,4/0.5,8/1"    ],
  300:   [ 0,     0,        77, "2/0.1,10/0.2"       ],
  305:   [ 0,     0,        97, "2/0.1,5/0.2,6,10/2" ],
  310:   [ 0,     0,        97, "2/0.1,5/0.2,6,10/2" ],
  315:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  320:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  325:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  330:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  335:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  340:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  345:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  350:   [ 0,     0,        97, "2/0.1,3/0.2,4,10/2" ],
  360:   [ 0,     0,        97, "0.2,1.8/0.1,2,6/2"  ], # 0.2-1.8
  370:   [ 0,     0,        97, "2/0.1,3/0.2,4,10/2" ],
  380:   [ 0,     0,        97, "2/0.1,3/0.2,4,10/2" ],
  390:   [ 0,     0,        97, "2/0.1,3/0.2,4,10/2" ],
  400:   [ 0,     0,        25, "2/0.1,3/0.2,4,10/2" ],
  420:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  440:   [ 0,     0,        97, "2/0.1,4/0.2,10/2"   ],
  460:   [ 0,     0,        97, "0.2,2/0.1,2.8/0.2,4,10/2" ], # 0.3-2.8
  480:   [ 0,     0,        97, "2/0.1,5/0.2,6,10/2" ],
  500:   [ 0,     0,        97, "2/0.1,5/0.2,6,10/2" ],
  520:   [ 0,     0,        97, "2/0.1,6/0.2,10/2"   ],
  540:   [ 0,     0,        97, "2/0.1,6/0.2,10/2"   ],
  560:   [ 0,     0,        97, "0.6/0.2,2/0.1,5/0.2,6,10/2,15" ],     # 0.6-5.0
  580:   [ 0,     0,        97, "10/0.2,11/0.5,15/1" ],
  600:   [ 0,     0,        97, "10/0.2,11/0.5,15/1" ],
}
