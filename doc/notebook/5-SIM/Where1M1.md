# WHERE1-M1 UWB measurement campaign

```python
>>> %matplotlib inline
```

```python
>>> from pylayers.simul.simulem import *
>>> from pylayers.antprop.rays import *
>>> from pylayers.antprop.channel import *
>>> from pylayers.antprop.signature import *
>>> from pylayers.measures.mesuwb import *
>>> import pylayers.util.pyutil as pyu
>>> import pylayers.signal.bsignal as bs
WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.
```

The deliverable describing the FP7 WHERE1 measurement campaign M1 can be found @
[WHERE1 D4.1 Measurements of location-dependent channel features](http://www.kn-s.dlr.de/where/documents/Deliverable41.pdf)

## Simulation Creation

A ray-tracing simulation is controlled via a configuration file which is stored in the `ini` directory of the project directory.
By default, there is a default configuration file named `default.ini`.

A first step consists in loading a layout associated with the simulation.
Here the `WHERE1.ini` layout is chosen along with the corresponding slabs and materials.

This layout corresponds to the office building where the first WHERE1 UWB measurement campaign has been conducted.


The layout method loads those in a member layout object **L** of the simulation object **S**.

If not already available, the layout associated graphs are built.

```python
>>> S = Simul()
>>> # loading a layout
... filestr = 'WHERE1'
>>> S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
>>> try:
...     S.L.dumpr()
>>> except:
...     S.L.build()
...     S.L.dumpw()
new file WHERE1.str
```

```python
>>> f,a=S.L.showG('s',figsize=(16,8),nodes=False)
```

```python
>>> S.L._showGi(en=17)
>>> plt.axis('off')
int0 :  (220, 58, 61)
int1 :  (216, 61)
 output  {(223, 61): 4.631245566616953e-05, (189, 61): 0.52059781654978099, (213, 61): 0.99999999999997502, (222, 61): 0.42854915043759029, (207, 61, 60): 0.015898410593098666, (213, 61, 59): 0.99999999999997502, (207, 61): 0.015898410593098666, (223, 61, 57): 4.631245566616953e-05, (205, 61): 0.074214608852340755, (189, 61, 70): 0.52059781654978099, (205, 61, 60): 0.074214608852340755}
(-31.123950000000001, 34.74295, 3.6289500000000001, 24.594305147086835)
```

The layout display is fully parameterized via the embedded **display** dictionnary member of the Layout object.

```python
>>> fig = plt.figure(figsize=(10,5))
>>> S.L.display['ednodes']=False
>>> S.L.display['nodes']=False
>>> S.L.display['title']='WHERE1 Project Office Measurement Site'
>>> fig,ax=S.L.showGs(fig=fig)
```

```python
>>> S.L.Gi.edges()[0]
((274, 55, 54), (286, 52, 55))
```

## Adding coordinates of transmiting and receiving points

Transmitters and receivers coordinates for the simulation are stored in **.ini** files.
Transmitter and Receiver are instances of the class `RadioNode` which offers different methods
for specifying nodes positions.
The stucture of this **.ini** file presented below.
The node Id is associated with the 3 coordinates $x,y,z$ separated by white spaces.

    [coordinates]
    1 = -12.2724 7.76319999993 1.2
    2 = -18.7747 15.1779999998 1.2
    3 = -4.14179999998 8.86029999983 1.2
    4 = -9.09139999998 15.1899000001 1.2

```python
>>> S.tx = RadioNode(_fileini='w2m1rx.ini',_fileant='defant.vsh3')
>>> S.rx = RadioNode(_fileini='w2m1tx.ini',_fileant='defant.vsh3')
```

The whole simulation setup can then be displayed using the **show** method of the Simulation object

```python
>>> fig = plt.figure(figsize=(10,5))
>>> fig,ax = S.show()
Warning : no furniture file loaded
```

Select Tx and Rx positions

```python
>>> map={1: 1,2: 2, 3: 3, 4: 5, 5: 6, 6: 7, 7: 8, 8: 9, 9: 10,
>>> 10: 11, 11: 12, 12: 13, 13: 14, 14: 15, 15: 16, 16: 17, 17: 18, 18: 19, 19: 20,
>>> 20: 21, 21: 22, 22: 23, 23: 24, 24: 25, 25: 26, 26: 27,
>>> 27: 28, 28: 29, 29: 30, 30: 32, 31: 33, 32: 34, 33: 35, 34: 36, 35: 37, 36: 38,
...       37: 39, 38: 40, 39: 41, 40: 42, 41: 43, 42: 44, 43: 45, 44: 46, 45: 47,
...       46: 48, 47: 49, 48: 50, 49: 51, 50: 52, 51: 53, 52: 54, 53: 55, 54: 56,
...       55: 57, 56: 58, 57: 59, 58: 60, 59: 61, 60: 62, 61: 63, 62: 64, 63: 65,
...       64: 66, 65: 67, 66: 68, 67: 69, 68: 70, 69: 71, 70: 72, 71: 73, 72: 74,
...       73: 75, 74: 76, 75: 77, 76: 78, 77: 79, 78: 80, 79: 81, 80: 82, 81: 83,
...       82: 84, 83: 85, 84: 89, 85: 90, 86: 91, 87: 92, 88: 93, 89: 94, 90: 95,
...       91: 96, 92: 97, 93: 98, 94: 99, 95: 100, 96: 101, 97: 103, 98: 104, 99:
...       105, 100: 106, 101: 107, 102: 108, 103: 109, 104: 110, 105: 111, 106:
...       113, 107: 114, 108: 116, 109: 117, 110: 119, 111: 120, 112: 122, 113:
...       123, 114: 124, 115: 125, 116: 126, 117: 127, 118: 128, 119: 129, 120:
...       133, 121: 134, 122: 136, 123: 137, 124: 138, 125: 139, 126: 140, 127:
...       141, 128: 142, 129: 143, 130: 144, 131: 145, 132: 146, 133: 147, 134:
...       162, 135: 163, 136: 164, 137: 165, 138: 166, 139: 167, 140: 168, 141:
...       169, 142: 170, 143: 171, 144: 172, 145: 173, 146: 174, 147: 175, 148:
...       176, 149: 177, 150: 179, 151: 180, 152: 181, 153: 182, 154: 183, 155:
...       184, 156: 185, 157: 186, 158: 188, 159: 189, 160: 199, 161: 200, 162:
...       201, 163: 202, 164: 203, 165: 204, 166: 205, 167: 206, 168: 207, 169:
...       208, 170: 209, 171: 210, 172: 211, 173: 212, 174: 213, 175: 214, 176:
...       215, 177: 216, 178: 217, 179: 218, 180: 219, 181: 220, 182: 221, 183:
...       222, 184: 223, 185: 227, 186: 228, 187: 229, 188: 230, 189: 231, 190:
...       232, 191: 233, 192: 234, 193: 235, 194: 236, 195: 237, 196: 238, 197:
...       239, 198: 240, 199: 241, 200: 242, 201: 243, 202: 244, 203: 245, 204:
...       246, 205: 247, 206: 248, 207: 249, 208: 250, 209: 251, 210: 252, 211:
...       253, 212: 258, 213: 259, 214: 266, 215: 267, 216: 268, 217: 269, 218:
...       270, 219: 271, 220: 272, 221: 273, 222: 274, 223: 275, 224: 276, 225:
...       277, 226: 278, 227: 279, 228: 297, 229: 298, 230: 299, 231: 300, 232:
...       301, 233: 302, 234: 303, 235: 304, 236: 305, 237: 306, 238: 307, 239:
...       308, 240: 309, 241: 310, 242: 311, 243: 312, 244: 313, 245: 314, 246:
...       315, 247: 316, 248: 317, 249: 318, 250: 319, 251: 320, 252: 321, 253:
...       322, 254: 323, 255: 324, 256: 325, 257: 326, 258: 327, 259: 328, 260:
...       329, 261: 330, 262: 332, 263: 333, 264: 334, 265: 335, 266: 336, 267:
...       337, 268: 338, 269: 339, 270: 340, 271: 341, 272: 342, 273: 343, 274:
...       344, 275: 345, 276: 346, 277: 347, 278: 348, 279: 349, 280: 350, 281:
...       351, 282: 352, 283: 353, 284: 354, 285: 355, 286: 356, 287: 360, 288:
...       361, 289: 362, 290: 363, 291: 364, 292: 365, 293: 366, 294: 367, 295:
...       368, 296: 369, 297: 370, 298: 371, 299: 372, 300: 373, 301: 374, 302:
...       375}
```

```python
>>> print 'number of Tx :',len(S.tx.points.keys())
>>> print 'number of rx :',len(S.rx.points.keys())
number of Tx : 2
number of rx : 350
```

Choose measurement points

```python
>>> # Chose used points here
... itx=10
>>> irx=2
>>> # check points
... tx= S.tx.points[itx]
>>> rx= S.rx.points[irx]
>>> M = UWBMeasure(map[itx])
>>> txm = M.tx
>>> rxm = M.rx[irx]
>>> print tx,txm
>>> print rx,rxm
>>> 
>>> if (tx[0] - txm[0] > 0.001) or (tx[1] - txm[1] > 0.001):
...     print 'Tx and Txm are not the same !'
>>> else :
...     print 'Txs OK'
>>> if (rx[0] - rxm[0] > 0.001) or (rx[1] - rxm[1] > 0.001):
...     print 'Rx and Rxm are not the same !'
>>> else :
...     print 'Rxs OK'
[-24.867   12.3097   1.2   ] [-24.867   12.3097   1.2   ]
[-18.7747  15.178    1.2   ] [-18.7747  15.178    1.2   ]
Txs OK
Rxs OK
```

```python
>>> fig =plt.figure(figsize=(16,8))
>>> fig,ax=S.L.showG('s',fig=fig)
>>> ax.plot(M.tx[0],M.tx[1],'or',label='tx')
>>> ax.plot(M.rx[irx][0],M.rx[irx][1],'ob',label='rx')
>>> ax.legend()
```

## Signatures, Rays and Radio Channel

A signature is a sequence of layout objects (points and segments) which are involved in a given optical ray joint the transmiter and the receiver.
The signatutre is calculated from a layout cycle to an other layout cycle. This means that is is required first to retrieve the cycle number from
point coordinates. This is done thanks to the **pt2cy**, point to cycle function.

```python
>>> ctx=S.L.pt2cy(tx)
>>> crx=S.L.pt2cy(rx)
>>> print 'tx point belongs to cycle ',ctx
>>> print 'rx point belongs to cycle ',crx
tx point belongs to cycle  7
rx point belongs to cycle  4
```

Then the signature between 2 given cycle can be calculated. This is done by instantiating a Signature object with a given layout and the 2 cycle number.

The representaion of a signature objet

```python
>>> Si = Signatures(S.L,ctx,crx)
>>> Si.run5(cutoff=3)
```

```python
>>> tx[2]=1.5
```

```python
>>> r2d = Si.rays(tx,rx)
>>> r3d = r2d.to3D(S.L)
```

```python
>>> fig = plt.figure(figsize=(10,10))
>>> r2d.show(L=S.L,fig=fig)
(<matplotlib.figure.Figure at 0x7f876d481ad0>,
 <matplotlib.axes.AxesSubplot at 0x7f876cd4b850>)
```

```python
>>> r3d.locbas(S.L)
>>> r3d.fillinter(S.L)
>>> r3d
Rays3D
----------
1 / 1 : [0]
2 / 6 : [1 2 3 4 5 6]
3 / 32 : [ 7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
 32 33 34 35 36 37 38]
4 / 137 : [ 39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56
  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74
  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92
  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110
 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128
 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146
 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164
 165 166 167 168 169 170 171 172 173 174 175]
5 / 227 : [176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193
 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211
 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229
 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247
 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265
 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283
 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301
 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319
 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337
 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355
 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373
 374 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391
 392 393 394 395 396 397 398 399 400 401 402]
6 / 198 : [403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420
 421 422 423 424 425 426 427 428 429 430 431 432 433 434 435 436 437 438
 439 440 441 442 443 444 445 446 447 448 449 450 451 452 453 454 455 456
 457 458 459 460 461 462 463 464 465 466 467 468 469 470 471 472 473 474
 475 476 477 478 479 480 481 482 483 484 485 486 487 488 489 490 491 492
 493 494 495 496 497 498 499 500 501 502 503 504 505 506 507 508 509 510
 511 512 513 514 515 516 517 518 519 520 521 522 523 524 525 526 527 528
 529 530 531 532 533 534 535 536 537 538 539 540 541 542 543 544 545 546
 547 548 549 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564
 565 566 567 568 569 570 571 572 573 574 575 576 577 578 579 580 581 582
 583 584 585 586 587 588 589 590 591 592 593 594 595 596 597 598 599 600]
7 / 30 : [601 602 603 604 605 606 607 608 609 610 611 612 613 614 615 616 617 618
 619 620 621 622 623 624 625 626 627 628 629 630]
8 / 4 : [631 632 633 634]
-----
ni : 3222
nl : 7079
```

```python
>>> S.freq()[0:10]
array([ 2.  ,  2.05,  2.1 ,  2.15,  2.2 ,  2.25,  2.3 ,  2.35,  2.4 ,  2.45])
```

```python
>>> Ct = r3d.eval(S.freq())
```

The `energy` method calculates the energy of each ray

```python
>>> Ett,Epp,Etp,Ept = Ct.energy()
```

```python
>>> plt.subplot(121)
>>> plt.plot(Ct.tauk,10*np.log10(Ett),'ob',label=r'$\theta\theta$')
>>> plt.plot(Ct.tauk,10*np.log10(Epp),'or',label=r'$\phi\phi$')
>>> plt.ylim(-160,-60)
>>> plt.xlabel('delay(ns)')
>>> plt.ylabel('Ray Energy (dB)')
>>> plt.legend()
>>> plt.subplot(122)
>>> plt.plot(Ct.tauk,10*np.log10(Ept),'og',label =r'$\phi\theta$')
>>> plt.plot(Ct.tauk,10*np.log10(Etp),'oc',label = r'$\theta\phi$')
>>> plt.ylim(-160,-60)
>>> plt.legend()
>>> plt.xlabel('delay(ns)')
```

## Apply waveform

```python
>>> Ct.freq = S.freq
>>> sco= Ct.prop2tran(a='theta',b='theta')
>>> sca= Ct.prop2tran(a=S.tx.A,b=S.rx.A)
```

```python
>>> wav = wvf.Waveform(typ='W1offset')
>>> #wav = wvf.Waveform({'type' : 'generic','band': 0.499,'fc': 4.493, 'fe': 100, 'thresh': 3, 'tw': 30})
... wav.show()
```

```python
>>> sco.isFriis
True
```

```python
>>> if sco.isFriis:
...     ciro = sco.applywavB(wav.sf)
>>> else:
...     ciro = sco.applywavB(wav.sfg)
>>> if sca.isFriis:
...     cira = sca.applywavB(wav.sf)
>>> else:
...      cira = sca.applywavB(wav.sfg)
```

```python
>>> ciro.plot(typ='v')
>>> f=plt.title(u'received waveform without antenna $\\theta\\theta$')
```

```python
>>> cira.plot(typ='v')
>>> f=plt.title('received waveform with antenna')
```

```python
>>> #dchan={i:'ch'+str(i) for i in range(1,5)}
... dchan={}
>>> dchan[1]='ch3'
>>> dchan[2]='ch4'
>>> dchan[3]='ch1'
>>> dchan[4]='ch2'
```

```python
>>> fig = plt.figure(figsize=(10,6))
>>> ax1 = fig.add_subplot(311,title="Measurements")
>>> cmd='M.tdd.' + str(dchan[irx]) + '.plot(ax=ax1)'
>>> eval(cmd)
>>> plt.title('WHERE1 measurement')
>>> #M.tdd.ch2.plot()
... # align for plotting
... #ciro.x=ciro.x-ciro.x[0]
... ax2 = fig.add_subplot(312,title="Simulation-with antenna",sharex=ax1, sharey=ax1)
>>> plt.xlim(20,70)
>>> plt.ylim(-95,-50)
>>> u = cira.plot(ax=ax2)
>>> plt.title('Simulation-with antenna - without noise')
>>> plt.tight_layout()
>>> #ax3 = fig.add_subplot(313,title="Simulation-without antenna",sharex=ax1, sharey=ax1)
... #ciro.plot()
```

```python
>>> r3d.info(0)
-------------------------
Informations of ray # 0
-------------------------

Index , type, slab      , th(rad), alpha     , gamma2    
    0 , B0  , -         , -      , -         , -         
    0 , T   , PARTITION ,    0.43,  (0.35+0j),  (0.88+0j)
    0 , B   , -         , -      , -         , -         

----------------------------------------
 Matrix of ray # 0 at f= 2.0
----------------------------------------
rotation matrix# type: B0
[[-0.99533359 -0.09649373]
 [ 0.09649373 -0.99533359]]
interaction # 0 type: T
[[ 0.05599940-0.55542654j  0.00000000-0.j        ]
 [ 0.00000000-0.j          0.03832677-0.60999405j]]
rotation matrix# [0] type: B
[[-0.99533359 -0.09649373]
 [ 0.09649373 -0.99533359]]
```

```python
>>> f,a=Ct.doadod(phi=(-180,180),cmap='copper')
```
