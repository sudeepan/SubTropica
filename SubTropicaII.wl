(* ::Package:: *)

BeginPackage["SubTropicaII`",{"SubTropica`","HyperIntica`","Parallel`Developer`"}];


Print["SubTropicaII Loaded!"]


(* ::Section:: *)
(*Preamble*)


(* ::Section::Closed:: *)
(*Description*)


STIIhyperIntMaple::usage=""
STsuggestOrderMaple::usage=""
STcgReductionMaple::usage=""


(* ::Section::Closed:: *)
(*Private*)


(*Begin["`Private`"];*)


(* ::Section:: *)
(*Tropical Subtractions II*)


(* ::Subsection::Closed:: *)
(*u & v-functions*)


(* Sort w-vectors by their complexity *)
STScoreWvecs[v_]:={Max[Abs/@v],Count[Abs/@v,Max[Abs/@v]],Total[Abs/@v]}


(* Produce the simplest u-variables for given rays and divergent faces *)
Clear[STProduceUs];
STProduceUs[xvars_,rays_,divFaces_,asAnAssociationQ_:False,normalizeQ_:True,QuietQ_:True]:=Module[
	{
	genMonomial,exponents,conditionLinear,compDiv,incompDiv,
	firstNonZero,
	facets,pairsFacets,coxCoords,us,nullVector,nullSpace,geometricQ,nullSpaceTrivialQ
	},
	facets=Cases[divFaces,{_}]//Flatten;
	pairsFacets=Cases[divFaces,{_,_}];
	(* Pairs of compatible facets *)
	compDiv=Association@@Table[f-> (Cases[pairsFacets,{f,_}|{_,f}]//Flatten//DeleteCases[#,f]&//DeleteDuplicates),{f,facets}];
	(* Pairs of incompatible facets *)
	incompDiv=Association@@Table[f->(Complement[facets,Join[{f},compDiv[f]]]),{f,facets}];
	us=
	Table[
		If[compDiv[f]==={},
			(*nullVector=-rays[[f]];*)
			(* A simplest w-vector *)
			firstNonZero=rays[[f]]//FirstPosition[#,v_/;Not[v===0],Heads->False]&//First;
			nullVector=-Table[KroneckerDelta[v,firstNonZero]rays[[f,v]],{v,1,rays[[f]]//Length}];
			(* The geometric property always holds if there are no compatible facets *)
			geometricQ=True;
			,
			(* Possible w-vectors *)
			nullSpace=rays[[compDiv[f]]]//NullSpace//SortBy[#,STScoreWvecs]&;
			If[nullSpace==={},
				nullSpace={Table[0,rays//First//Length]};
				nullSpaceTrivialQ=True;
				,
				nullSpaceTrivialQ=False;
			];
			If[Not[QuietQ],
				Print["The space orthogonal to facet "<>ToString[f]<>" is spanned by: " <>ToString[nullSpace]];
			];
			(* Remove w-vectors orthogonal to facet *)
			If[Not[nullSpaceTrivialQ],
				nullSpace=nullSpace//DeleteCases[#,v_/;v . rays[[f]]===0]&;
			];
			nullVector=nullSpace//First; 
			(* Does the geometric property hold? *)
			geometricQ=Not[nullSpaceTrivialQ] && Not[nullSpace==={}];
		];
		(* Normalize the vector so that w.facet = - 1 *)
		If[geometricQ&&normalizeQ,
			nullVector=nullVector/(nullVector . rays[[f]])
		];
		(* If necessary flip the sign *)
		If[nullVector . rays[[f]]>0&&geometricQ,nullVector=-nullVector];
		(* If geometric property holds, construct v-function *)
		If[geometricQ,
			1-1/(1+Times@@(xvars^nullVector))//Factor
			,
			"NotFound"
		]
		,
		{f,facets[[;;]]}
	];
	(* If asked to, return an association *)
	If[asAnAssociationQ,
		Association@@Table[facets[[i]]-> us[[i]],{i,1,facets//Length}],
		us
	]
];


(* Produce all u-variables *)
Clear[STProduceAllUs];STProduceAllUs[xvars_,rays_,divFaces_,normalizeQ_:True]:=Module[{genMonomial,exponents,conditionLinear,compDiv,incompDiv,facets,pairsFacets,nullVectors,nullVector},
facets=Cases[divFaces,{_}]//Flatten;
pairsFacets=Cases[divFaces,{_,_}];
compDiv=Association@@Table[f-> (Cases[pairsFacets,{f,_}|{_,f}]//Flatten//DeleteCases[#,f]&//DeleteDuplicates),{f,facets}];
incompDiv=Association@@Table[f->(Complement[facets,Join[{f},compDiv[f]]]),{f,facets}];
Table[
nullVectors=If[compDiv[f]==={},
-rays[[{f}]],
rays[[compDiv[f]]]//NullSpace//SortBy[#,STScoreWvecs]&//DeleteCases[#,v_/;v . rays[[f]]===0]&
];
If[nullVectors==={},"NotFound",
Table[
nullVector=nV;
If[Not[nullVector . rays[[f]]===0]&&Not[nullVector===First[{}]]&&normalizeQ,
nullVector=nullVector/(nullVector . rays[[f]])
];
If[nullVector . rays[[f]]>0&&Not[nullVector===First[{}]],nullVector=-nullVector];
If[nullVector . rays[[f]]===0||nullVector===First[{}],"NotFound",1-1/(1+Times@@(xvars^nullVector))//Factor],{nV,nullVectors}]//DeleteDuplicates
]
,{f,facets[[;;]]}]/.{"NotFound"}-> "NotFound"
]


(* ::Subsection::Closed:: *)
(*Geometric Property Analysis*)


(* Returns divergent rays, trop values, faces, u-variables, trData *)
Clear[STIIPreAnalysis];
STIIPreAnalysis[integrand_,vars_,coeffs_,facesQ_:True,regulators_:{eps}]:=Module[
	{trData,regularizeFace,monspols,polsInt,divFacets,faces,
	divFaces,us,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,newIntegrand,result,A,B,trint,trValues,Wvects,
	wherePowerDiv,ordersOfDivergences	
	}
	,
		
	(* Computing Tropical Data *)
	Print["Computing tropical data"];
	trData=integrand//STIntegrandTropicalData[#,vars,coeffs]&;
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	If[facesQ,faces=trData//STGetFaces[#,divFacets]&,faces="NotComputed"];
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	ordersOfDivergences=Association@@Table[facet-> (trValues[[facet]]/.eps->0),{facet,divFacets}];
	
	If[MemberQ[Values[ordersOfDivergences],v_/;v>0],
		wherePowerDiv=divFacets[[ Position[Values[ordersOfDivergences],v_/;v>0]//Flatten ]];
		Print["Power divergences on rays "<>ToString[wherePowerDiv]];
	];
	
	(*u-variables*)
	us=Table[{},{f,trData["rays"]}];
	us[[divFacets]]=STProduceUs[vars,trData["rays"],divFaces]//Factor;
	Association@@{
		"trops"-> trValues[[divFacets]], 
		"rays" -> divFacets,
		"faces" -> divFaces, 
		"us" -> us[[divFacets]], 
	"trData"-> trData
	}
];


(* ::Subsection::Closed:: *)
(*Subtraction Formula [new version; NP updated]*)


(* Compute subtraction formula *) 
(* It is possible to pass tropical data and u-variables if these have already been found *)
Clear[STIISubtractionFormula];
Options[STIISubtractionFormula]={
	"Gauges"->"Automatic"
};
STIISubtractionFormula[integrand_,vars_,coeffs_,divUs_:{},regulators_:{eps},trDataGiven_:{}]:=Module[
	{trData,regularizeFace,monspols,polsInt,divFacets,faces,divFaces,
	us,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,
	newIntegrand,result,A,B,trint,trValues,Wvects},
	Print["Computing tropical data"];
	trData=If[trDataGiven==={},integrand//STIntegrandTropicalData[#,vars,coeffs]&,trDataGiven];
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	If[Count[trValues/.Alternatives@@regulators->0,v_/;v> 0]>0,
		Print["Power divergence on "<>ToString[ Position[trValues/.Alternatives@@regulators->0,v_/;v> 0]//Flatten];]
	];
	
	monspols=STtoMonPols[integrand,vars];
	polsInt=Times@@(STtoMonPols[integrand,vars][[2,1]]);
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	(* If no divergences, return integrand in subtraction formula format.
	   Divide by Times@@vars to convert from dx/x to flat dx convention,
	   matching what regularizeFace does for non-trivial faces (line ~2507). *)
	If[divFacets==={},
		Print[" Integrand is already locally finite! " ];
		Return[{{{integrand/(Times@@vars)},vars}}];
	];
	
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	
	(*u-variables*)
	us=Table[{},{f,trData["rays"]}];
	If[Head[divUs]===Association,
	Table[us[[dF]]=divUs[dF],{dF,divFacets}]
	,
	Print["Computing u-variables "];
	(* STProduceUs[xvars_,rays_,divFaces_,asAnAssociationQ_:False,normalizeQ_:True,QuietQ_:True] *)
	us[[divFacets]]=If[divUs==={},
	STProduceAllUs[vars,trData["rays"],divFaces,True]//Factor//Print;
	STProduceUs[vars,trData["rays"],divFaces,False,True,True]//Factor
	
	,divUs//Factor];
	
	];
	If[Count[us,"NotFound!"]>0,
		Print[" The geometric property is not satisfied! "];
		Return[us];
	,
	Print["Proceeding with u-variables: " , us ]
	];
	
	omus=1-us//Factor;
	
	(* Recover the W-vectors from the u-vars *)
	Wvects= Table[
	A = (1-us[[dF]])//Factor//Numerator;
	B =( (1-us[[dF]])//Factor//Denominator)-(A//Factor)//Factor;
	dF-> CoefficientRules[B,vars][[1,1]]-CoefficientRules[A,vars][[1,1]]
	,{dF,divFacets}];
	Wvects=Association@@Wvects;
	(*J space*)
	(*Js=Table[
	If[fc==={1,19},
	Print["special"];
	fc-> {2,3,4,6,7,8,9},
	fc-> FirstCase[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
	]
	,{fc,divFaces[[2;;]]}
	];
	Js=Join[{divFaces[[1]]-> Range[vars//Length]},Js];
	Js=Association@@Js;*)
	Js=Switch[OptionValue["Gauges"],
		"Automatic",
		Js=Table[
			fc-> FirstCase[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
			,
			{fc,divFaces[[2;;]]}
		];
		Js=Join[{divFaces[[1]]-> Range[vars//Length]},Js];
		Js=Association@@Js,
		
		"LRGauge",
		STIIFindLRGauges[integrand,xvars,coeffs]
	];
	Js//Print;
	
	(*Monomial*)
	extraMons=Association@@Table[
	fc->(1/(Times@@(vars^(-trValues[[fc]] . ( Wvects/@fc))))/.Alternatives@@vars[[Complement[Range[vars//Length],Js[fc]]]]->1)
	,{fc,divFaces}
	];
	
	(*Volume*)
	Print["Computing volumes"];
	volume=Table[
	fc->(1/Times@@(-trValues[[fc]]))(Abs[Det[Join[trData["rays"][[fc]],Table[KroneckerDelta[i,Js[fc][[j]]],{j,1,Length[Js[fc]]} ,{i,1,vars//Length}  ]]]])
	,{fc,divFaces[[2;;]]}
	];
	volume=Join[{divFaces[[1]]->1},volume];
	volume=Association@@volume;
	
	(*jacs = 1-u *)
	Print["Computing jacobians"];
	jacFace=Table[
	face->Table[Factor[(1-us[[j]])],{j,face}]
	,{face,divFaces}];
	jacFace = Association@@jacFace;
	
	(*We define a function which regularize a face*)
	regularizeFace[face_]:=Module[{compFacess,newIntegrandd},
	compFaces=Cases[divFaces,f_/;SubsetQ[f,face]];
	compFaces=Complement[#,face]&/@compFaces;
	(*We compute the face integrand to be regularized*)
	newIntegrand=If[face==={}, integrand,
	Times[
	STrestrictIntegrand[integrand,vars,trData["rays"][[face]]],
	extraMons[face]
	](*//STFactor//STSimplify*)
	];
	(*Print["computed new"];*)
	newIntegrand=newIntegrand;
	Monitor[
	result=volume[face] Table[
	If[subface==={},
	newIntegrand,
	Times[-a[subface],
	Times@@((jacFace[subface])^(1-trValues[[subface]]))(*why only 1+Trop, rather than n + Trop. Why set x[face]=0?*),
	STrestrictIntegrand[newIntegrand,vars,trData["rays"][[subface]]]
	]
	]
	,{subface,compFaces}];,
	{subface}];
	(*Print["computed result"];*)
	
	result=result(*//STFactor//STSimplify*);
	result=result/.a[subface_]:> (-1)^(Length[subface]+1);
	result=result/.Alternatives@@vars[[ Complement[Range[vars//Length],Js[face]] ]]->1;
	(*Print["computed gauge fix"];*)
	
	{result/(Times@@(vars[[Js[face]]])),vars[[Js[face]]]}
	];
	Print["Computing subtraction formula"];
	Monitor[Table[regularizeFace[divFaces[[face]] ],{face,1,divFaces//Length}]//Return,{face,"/",divFaces//Length}]
];

(* Alias requested by Giulio: STTropicalSubtraction mirrors STSubtractionFormula. *)
STTropicalSubtraction = STSubtractionFormula;


(* ::Subsection::Closed:: *)
(*Enumerate LR Gauge&Order*)


(* Enumerate Linearly Reducible Gauges *) 
(* The intended usage is by additionally passing a choice of v-functions (association <|ray->v_ray|>). Finally, regulators and trData can be optionally passed. *)
Clear[STIIEnumerateFaceLROrders];
STIIEnumerateFaceLROrders[integrand_,vars_,coeffs_,divUs_:{},regulators_:{eps},trDataGiven_:{}]:=Module[
	{trData,regularizeFace,monspols,polsInt,divFacets,faces,divFaces,
	us,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,
	newIntegrand,result,A,B,trint,trValues,Wvects,
	(* LR Gauges and Orders*)
	Iren,commonLRorders,LRorders,gauges,
	potentialLRorders,faceLRorders,
	polys,exps,LROrdersCompatibleWithGauges},
	
	Print["Computing tropical data"];
	trData=If[trDataGiven==={},integrand//STIntegrandTropicalData[#,vars,coeffs]&,trDataGiven];
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	
	If[Count[trValues/.Alternatives@@regulators->0,v_/;v> 0]>0,
		Print["Power divergence on "<>ToString[ Position[trValues/.Alternatives@@regulators->0,v_/;v> 0]//Flatten];]
	];
	
	monspols=STtoMonPols[integrand,vars];
	polsInt=Times@@(STtoMonPols[integrand,vars][[2,1]]);
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	(* If no divergences, return integrand in subtraction formula format.
	   Divide by Times@@vars to convert from dx/x to flat dx convention,
	   matching what regularizeFace does for non-trivial faces (line ~2507). *)
	If[divFacets==={},
		Print[" Note: Integrand is already locally finite! " ];
		{{},{},{}}//Return
	];
	
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	
	(*u-variables*)
	us=Table[{},{f,trData["rays"]}];
	If[Length[divUs]>0,
		Table[dF=divFacets[[j]];us[[dF]]=divUs[[j]],{j,1,divFacets//Length}]
	,
	Print["Computing u-variables "];
	(* STProduceUs[xvars_,rays_,divFaces_,asAnAssociationQ_:False,normalizeQ_:True,QuietQ_:True] *)
	us[[divFacets]]=If[divUs==={},
	(*STProduceAllUs[vars,trData["rays"],divFaces,True]//Factor//Print;*)
	STProduceUs[vars,trData["rays"],divFaces,False,True,True]//Factor
	,divUs//Factor];
	
	];
	If[Count[us,"NotFound!"]>0,
		Print[" The geometric property is not satisfied! "];
		Return[us];
	,
	Print["Proceeding with u-variables: " , us ]
	];
	
	omus=1-us//Factor;
	
	(* Recover the W-vectors from the u-vars *)
	Wvects= Table[
	A = (1-us[[dF]])//Factor//Numerator;
	B =( (1-us[[dF]])//Factor//Denominator)-(A//Factor)//Factor;
	dF-> CoefficientRules[B,vars][[1,1]]-CoefficientRules[A,vars][[1,1]]
	,{dF,divFacets}];
	Wvects=Association@@Wvects;
	(*J space --- dummy, to be compatible with older code *)
	Js=Table[
		fc-> FirstCase[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
		,
		{fc,divFaces[[2;;]]}
	];
	Js=Join[{divFaces[[1]]-> Range[vars//Length]},Js];
	Js=Association@@Js;
	
	(*Monomial*)
	extraMons=Association@@Table[
	fc->(1/(Times@@(vars^(-trValues[[fc]] . ( Wvects/@fc))))/.Alternatives@@vars[[Complement[Range[vars//Length],Js[fc]]]]->1)
	,{fc,divFaces}
	];
	
	(*Volume*)
	Print["Computing volumes"];
	volume=Table[
	fc->(1/Times@@(-trValues[[fc]]))(Abs[Det[Join[trData["rays"][[fc]],Table[KroneckerDelta[i,Js[fc][[j]]],{j,1,Length[Js[fc]]} ,{i,1,vars//Length}  ]]]])
	,{fc,divFaces[[2;;]]}
	];
	volume=Join[{divFaces[[1]]->1},volume];
	volume=Association@@volume;
	
	(*jacs = 1-u *)
	Print["Computing jacobians"];
	jacFace=Table[
	face->Table[Factor[(1-us[[j]])],{j,face}]
	,{face,divFaces}];
	jacFace = Association@@jacFace;
	
	
	(* ---------- Renormalization Map ------------ *)
	(*We define a function which regularize a face*)
	regularizeFace[face_]:=Module[{compFacess,newIntegrandd},
	compFaces=Cases[divFaces,f_/;SubsetQ[f,face]];
	compFaces=Complement[#,face]&/@compFaces;
	(*We compute the face integrand to be regularized*)
	newIntegrand=If[face==={}, integrand,
	Times[
	STrestrictIntegrand[integrand,vars,trData["rays"][[face]]],
	extraMons[face]
	](*//STFactor//STSimplify*)
	];
	(*Print["computed new"];*)
	newIntegrand=newIntegrand;
	Monitor[
	result=volume[face] Table[
	If[subface==={},
	newIntegrand,
	Times[-a[subface],
	Times@@((jacFace[subface])^(1-trValues[[subface]]))(*why only 1+Trop, rather than n + Trop. Why set x[face]=0?*),
	STrestrictIntegrand[newIntegrand,vars,trData["rays"][[subface]]]
	]
	]
	,{subface,compFaces}];,
	{subface}];
	(*Print["computed result"];*)
	
	result=result(*//STFactor//STSimplify*);
	result=result/.a[subface_]:> (-1)^(Length[subface]+1);
	result=result/.Alternatives@@vars[[ Complement[Range[vars//Length],Js[face]] ]]->1;
	(*Print["computed gauge fix"];*)
	
	{result/(Times@@(vars[[Js[face]]])),vars[[Js[face]]]}
	];
	(* ---------- Renormalization Map ------------ *)
	
	(* --------- Possible Gauges --------- *)
	gauges=STIIFindGauges[integrand,vars,coeffs,trData,divFaces];
	
	(* ---------- Recursively Enumerate LR Gauges and orders ------------*)
	Iren=regularizeFace[divFaces[[1]]][[1]];
	(*{divFaces,Iren,LRorders,gauges//Values}//Return;*)
	
	(* Store LR orders for each un-integrated counter-term *)
	(* Better way: first find only common orders, then only orders ending with any of the gauges *)
	(*commonLRorders=Table[
		{polys,exps}=Iren[[i]]//STtoCoeffMonPols[#,vars]&//Last;
		polys=Table[If[Count[{exps[[i]]},eps,\[Infinity]]>0||exps[[i]]<0,polys[[i]],{}],{i,1,exps//Length}]//Flatten;
		polys
		,
		{i,1,divFaces//Length}
	]//STIIMapleFindLROrders[#,vars,coeffs]&;
	commonLRorders=commonLRorders//SortBy[#,#[[2]]&]&;*)
	LRorders=Table[
		{polys,exps}=Iren[[i]]//STtoCoeffMonPols[#,vars]&//Last;
		polys=Table[If[Count[{exps[[i]]},eps,\[Infinity]]>0||exps[[i]]<0,polys[[i]],{}],{i,1,exps//Length}]//Flatten;
		(*(polys//dummySTIIMapleFindAllLROrdersEndingWith[#,vars,coeffs,Values[gauges][[i]]]&)//Print;*)
		LROrdersCompatibleWithGauges=(polys//STIIMapleFindAllLROrdersEndingWith[#,vars,coeffs,Complement[vars,#]&/@(Values[gauges][[i]])]&);
		Print[ToString[i]<>"/"<>ToString[(divFaces//Length)-1]<>"+1"];
		LROrdersCompatibleWithGauges//Length//Print;
		LROrdersCompatibleWithGauges=Join[LROrdersCompatibleWithGauges]//SortBy[#,#[[2]]&]&;
		LROrdersCompatibleWithGauges[[;;,1]]
		If[Length[LROrdersCompatibleWithGauges]===0,Print["WARNING: no LR order compatible with gauge found for face: " <>ToString[i]]];
		
		,
		{i,2,divFaces//Length}
	];
	
	(* Prepend LR orders for integrand *)
	Print["last one: ..."];
	LRorders=Join[{commonLRorders},LRorders];
	(*LRorders=Join[{(integrand//STtoCoeffMonPols[#,vars]&//Last//First//STIIMapleFindLROrders[{#},vars,coeffs]&)[[;;,1]]},LRorders];*)
	(*Return[LRorders];*)
	backup5={divFaces,Iren,LRorders,gauges//Values};
	
	(**)
	Join[{commonLRorders[[;;,1]]},Table[
	potentialLRorders=Intersection@@(LRorders[[Position[divFaces,Alternatives@@Cases[divFaces,v_/;SubsetQ[v,divFaces[[face]] ]],1]//Flatten]])[[;;,;;,;;(Length[vars]-Length[divFaces[[face]]])]];
	faceLRorders=potentialLRorders//Cases[#,v_/;Or@@(Table[Sort[v]===Sort[ga],{ga,Values[gauges][[face]]}])]&
	,{face,2,divFaces//Length}]]
	
];



(* ::Subsection::Closed:: *)
(*Enumerate LR Gauge&Order 2*)


(*(* Enumerate Linearly Reducible Gauges *) 
(* The intended usage is by additionally passing a choice of v-functions (association <|ray->v_ray|>). Finally, regulators and trData can be optionally passed. *)
Clear[STIIEnumerateFaceLROrders2];
STIIEnumerateFaceLROrders2[integrand_,vars_,coeffs_,trDataGiven_:{}]:=Module[
	{trData,regularizeFace,monspols,polsInt,divFacets,faces,divFaces,
	us,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,
	newIntegrand,result,A,B,trint,trValues,Wvects,
	(* LR Gauges and Orders*)
	Iren,commonLRorders,LRorders,gauges,
	potentialLRorders,faceLRorders,
	polys,exps,LROrdersCompatibleWithGauges},
	
	Print["Computing tropical data"];
	trData=If[trDataGiven==={},integrand//STIntegrandTropicalData[#,vars,coeffs]&,trDataGiven];
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	
	If[Count[trValues/.Alternatives@@regulators->0,v_/;v> 0]>0,
		Print["Power divergence on "<>ToString[ Position[trValues/.Alternatives@@regulators->0,v_/;v> 0]//Flatten];]
	];
	
	monspols=STtoMonPols[integrand,vars];
	polsInt=Times@@(STtoMonPols[integrand,vars][[2,1]]);
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	(* If no divergences, return integrand in subtraction formula format.
	   Divide by Times@@vars to convert from dx/x to flat dx convention,
	   matching what regularizeFace does for non-trivial faces (line ~2507). *)
	If[divFacets==={},
		Print[" Note: Integrand is already locally finite! " ];
		{{},{},{}}//Return
	];
	
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	
	(*dummy u-variables*)
	us=Table[{},{f,trData["rays"]}];
	If[Length[divUs]>0,
		Table[dF=divFacets[[j]];us[[dF]]=0,{j,1,divFacets//Length}]
	,
	
	(* STProduceUs[xvars_,rays_,divFaces_,asAnAssociationQ_:False,normalizeQ_:True,QuietQ_:True] *)
	us[[divFacets]]=If[divUs==={},
	(*STProduceAllUs[vars,trData["rays"],divFaces,True]//Factor//Print;*)
	1
	,0];
	
	];
	If[Count[us,"NotFound!"]>0,
		Print[" The geometric property is not satisfied! "];
		Return[us];
	,
	Print["Proceeding with u-variables: " , us ]
	];
	
	omus=1-us//Factor;
	
	(* Recover the W-vectors from the u-vars *)
	Wvects= Table[
	(*A = (1-us[[dF]])//Factor//Numerator;*)
	(*B =( (1-us[[dF]])//Factor//Denominator)-(A//Factor)//Factor;*)
	dF-> 0
	,{dF,divFacets}];
	Wvects=Association@@Wvects;
	(*J space --- dummy, to be compatible with older code *)
	Js=Table[
		fc-> FirstCase[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
		,
		{fc,divFaces[[2;;]]}
	];
	Js=Join[{divFaces[[1]]-> Range[vars//Length]},Js];
	Js=Association@@Js;
	
	(*Monomial*)
	extraMons=Association@@Table[
	fc->(1/(Times@@(vars^(-trValues[[fc]] . ( Wvects/@fc))))/.Alternatives@@vars[[Complement[Range[vars//Length],Js[fc]]]]->1)
	,{fc,divFaces}
	];
	
	(*Volume*)
	Print["Computing volumes"];
	volume=Table[
		fc->1
		,
		{fc,divFaces[[2;;]]}
	];
	volume=Join[{divFaces[[1]]->1},volume];
	volume=Association@@volume;
	
	(*jacs = 1-u *)
	Print["Computing jacobians"];
	jacFace=Table[
		face->1
		,
		{face,divFaces}
	];
	jacFace = Association@@jacFace;
	
	(* ---------- Renormalization Map ------------ *)
	(*We define a function which regularize a face*)
	regularizeFace[face_]:=Module[{compFacess,newIntegrandd},
	compFaces=Cases[divFaces,f_/;SubsetQ[f,face]];
	compFaces=Complement[#,face]&/@compFaces;
	(*We compute the face integrand to be regularized*)
	newIntegrand=If[face==={}, integrand,
		Times[
		STrestrictIntegrand[integrand,vars,trData["rays"][[face]]]
		,
		extraMons[face]
		]
	];
	(*Print["computed new"];*)
	newIntegrand=newIntegrand;
	Monitor[
	result=volume[face] Table[
	If[subface==={},
	newIntegrand,
	Times[-a[subface],
	Times@@((jacFace[subface])^(1-trValues[[subface]]))(*why only 1+Trop, rather than n + Trop. Why set x[face]=0?*),
	STrestrictIntegrand[newIntegrand,vars,trData["rays"][[subface]]]
	]
	]
	,{subface,compFaces}];,
	{subface}];
	(*Print["computed result"];*)
	
	result=result(*//STFactor//STSimplify*);
	result=result/.a[subface_]:> (-1)^(Length[subface]+1);
	result=result/.Alternatives@@vars[[ Complement[Range[vars//Length],Js[face]] ]]->1;
	(*Print["computed gauge fix"];*)
	
	{result/(Times@@(vars[[Js[face]]])),vars[[Js[face]]]}
	];
	(* ---------- Renormalization Map ------------ *)
	
	(* --------- Possible Gauges --------- *)
	gauges=STIIFindGauges[integrand,vars,coeffs,trData,divFaces];
	Table[
	
	(* ---------- Recursively Enumerate LR Gauges and orders ------------*)
	Iren=regularizeFace[divFaces[[1]]][[1]];
	(*{divFaces,Iren,LRorders,gauges//Values}//Return;*)
	
	(* Store LR orders for each un-integrated counter-term *)
	(* Better way: first find only common orders, then only orders ending with any of the gauges *)
	(*commonLRorders=Table[
		{polys,exps}=Iren[[i]]//STtoCoeffMonPols[#,vars]&//Last;
		polys=Table[If[Count[{exps[[i]]},eps,\[Infinity]]>0||exps[[i]]<0,polys[[i]],{}],{i,1,exps//Length}]//Flatten;
		polys
		,
		{i,1,divFaces//Length}
	]//STIIMapleFindLROrders[#,vars,coeffs]&;
	commonLRorders=commonLRorders//SortBy[#,#[[2]]&]&;*)
	
	
	LRorders=Table[
		{polys,exps}=Iren[[i]]//STtoCoeffMonPols[#,vars]&//Last;
		polys=Table[If[Count[{exps[[i]]},eps,\[Infinity]]>0||exps[[i]]<0,polys[[i]],{}],{i,1,exps//Length}]//Flatten;
		(*(polys//dummySTIIMapleFindAllLROrdersEndingWith[#,vars,coeffs,Values[gauges][[i]]]&)//Print;*)
		LROrdersCompatibleWithGauges=(polys//STIIMapleFindAllLROrdersEndingWith[#,vars,coeffs,Complement[vars,#]&/@(Values[gauges][[i]])]&);
		Print[ToString[i]<>"/"<>ToString[(divFaces//Length)-1]<>"+1"];
		LROrdersCompatibleWithGauges//Length//Print;
		LROrdersCompatibleWithGauges=Join[LROrdersCompatibleWithGauges]//SortBy[#,#[[2]]&]&;
		If[Length[LROrdersCompatibleWithGauges]===0,Print["WARNING: no LR order compatible with gauge found for face: " <>ToString[i]]];
		LROrdersCompatibleWithGauges[[;;,1]];
		,
		{i,2,divFaces//Length}
	];
	
	(* Prepend LR orders for integrand *)
	Print["last one: ..."];
	LRorders=Join[{commonLRorders},LRorders];
	(*LRorders=Join[{(integrand//STtoCoeffMonPols[#,vars]&//Last//First//STIIMapleFindLROrders[{#},vars,coeffs]&)[[;;,1]]},LRorders];*)
	(*Return[LRorders];*)
	backup5={divFaces,Iren,LRorders,gauges//Values};
	
	(**)
	Join[{commonLRorders[[;;,1]]},Table[
	potentialLRorders=Intersection@@(LRorders[[Position[divFaces,Alternatives@@Cases[divFaces,v_/;SubsetQ[v,divFaces[[face]] ]],1]//Flatten]])[[;;,;;,;;(Length[vars]-Length[divFaces[[face]]])]];
	faceLRorders=potentialLRorders//Cases[#,v_/;Or@@(Table[Sort[v]===Sort[ga],{ga,Values[gauges][[face]]}])]&
	,{face,2,divFaces//Length}]]
	
];


*)


(* ::Subsection::Closed:: *)
(*Enumerate LR Gauge&Order 3*)


(*(* Enumerate Linearly Reducible Gauges *) 
(* The intended usage is by additionally passing a choice of v-functions (association <|ray->v_ray|>). Finally, regulators and trData can be optionally passed. *)
Clear[STIIEnumerateFaceLROrders3];
STIIEnumerateFaceLROrders3[integrand_,vars_,coeffs_,trDataGiven_:{}]:=Module[
	{trData,regularizeFace,monspols,polsInt,divFacets,faces,divFaces,
	us,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,
	newIntegrand,result,A,B,trint,trValues,Wvects,
	(* LR Gauges and Orders*)
	Iren,commonLRorders,LRorders,gauges,
	potentialLRorders,faceLRorders,
	polys,exps,LROrdersCompatibleWithGauges},
	
	Print["Computing tropical data"];
	trData=If[trDataGiven==={},integrand//STIntegrandTropicalData[#,vars,coeffs]&,trDataGiven];
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	
	If[Count[trValues/.Alternatives@@regulators->0,v_/;v> 0]>0,
		Print["Power divergence on "<>ToString[ Position[trValues/.Alternatives@@regulators->0,v_/;v> 0]//Flatten];]
	];
	
	monspols=STtoMonPols[integrand,vars];
	polsInt=Times@@(STtoMonPols[integrand,vars][[2,1]]);
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	(* If no divergences, return integrand in subtraction formula format.
	   Divide by Times@@vars to convert from dx/x to flat dx convention,
	   matching what regularizeFace does for non-trivial faces (line ~2507). *)
	If[divFacets==={},
		Print[" Note: Integrand is already locally finite! " ];
		{{},{},{}}//Return
	];
	
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	
	(*dummy u-variables*)
	us=Table[{},{f,trData["rays"]}];
	If[Length[divUs]>0,
		Table[dF=divFacets[[j]];us[[dF]]=0,{j,1,divFacets//Length}]
	,
	
	(* STProduceUs[xvars_,rays_,divFaces_,asAnAssociationQ_:False,normalizeQ_:True,QuietQ_:True] *)
	us[[divFacets]]=If[divUs==={},
	(*STProduceAllUs[vars,trData["rays"],divFaces,True]//Factor//Print;*)
	1
	,0];
	
	];
	If[Count[us,"NotFound!"]>0,
		Print[" The geometric property is not satisfied! "];
		Return[us];
	,
	Print["Proceeding with u-variables: " , us ]
	];
	
	omus=1-us//Factor;
	
	(* Recover the W-vectors from the u-vars *)
	Wvects= Table[
	(*A = (1-us[[dF]])//Factor//Numerator;*)
	(*B =( (1-us[[dF]])//Factor//Denominator)-(A//Factor)//Factor;*)
	dF-> 0
	,{dF,divFacets}];
	Wvects=Association@@Wvects;
	(*J space --- dummy, to be compatible with older code *)
	Js=Table[
		fc-> FirstCase[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
		,
		{fc,divFaces[[2;;]]}
	];
	Js=Join[{divFaces[[1]]-> Range[vars//Length]},Js];
	Js=Association@@Js;
	
	(*Monomial*)
	extraMons=Association@@Table[
	fc->(1/(Times@@(vars^(-trValues[[fc]] . ( Wvects/@fc))))/.Alternatives@@vars[[Complement[Range[vars//Length],Js[fc]]]]->1)
	,{fc,divFaces}
	];
	
	(*Volume*)
	Print["Computing volumes"];
	volume=Table[
		fc->1
		,
		{fc,divFaces[[2;;]]}
	];
	volume=Join[{divFaces[[1]]->1},volume];
	volume=Association@@volume;
	
	(*jacs = 1-u *)
	Print["Computing jacobians"];
	jacFace=Table[
		face->1
		,
		{face,divFaces}
	];
	jacFace = Association@@jacFace;
	
	(* ---------- Renormalization Map ------------ *)
	(*We define a function which regularize a face*)
	regularizeFace[face_]:=Module[{compFacess,newIntegrandd},
	compFaces=Cases[divFaces,f_/;SubsetQ[f,face]];
	compFaces=Complement[#,face]&/@compFaces;
	(*We compute the face integrand to be regularized*)
	newIntegrand=If[face==={}, integrand,
		Times[
		STrestrictIntegrand[integrand,vars,trData["rays"][[face]]]
		,
		extraMons[face]
		]
	];
	(*Print["computed new"];*)
	newIntegrand=newIntegrand;
	Monitor[
	result=volume[face] Table[
	If[subface==={},
	newIntegrand,
	Times[-a[subface],
	Times@@((jacFace[subface])^(1-trValues[[subface]]))(*why only 1+Trop, rather than n + Trop. Why set x[face]=0?*),
	STrestrictIntegrand[newIntegrand,vars,trData["rays"][[subface]]]
	]
	]
	,{subface,compFaces}];,
	{subface}];
	(*Print["computed result"];*)
	
	result=result(*//STFactor//STSimplify*);
	result=result/.a[subface_]:> (-1)^(Length[subface]+1);
	result=result/.Alternatives@@vars[[ Complement[Range[vars//Length],Js[face]] ]]->1;
	(*Print["computed gauge fix"];*)
	
	{result/(Times@@(vars[[Js[face]]])),vars[[Js[face]]]}
	];
	(* ---------- Renormalization Map ------------ *)
	
	(* --------- Possible Gauges --------- *)
	gauges=STIIFindGauges[integrand,vars,coeffs,trData,divFaces];
	Table[
	
	(* ---------- Recursively Enumerate LR Gauges and orders ------------*)
	Iren=regularizeFace[divFaces[[1]]][[1]];
	(*{divFaces,Iren,LRorders,gauges//Values}//Return;*)
	
	(* Store LR orders for each un-integrated counter-term *)
	(* Better way: first find only common orders, then only orders ending with any of the gauges *)
	(*commonLRorders=Table[
		{polys,exps}=Iren[[i]]//STtoCoeffMonPols[#,vars]&//Last;
		polys=Table[If[Count[{exps[[i]]},eps,\[Infinity]]>0||exps[[i]]<0,polys[[i]],{}],{i,1,exps//Length}]//Flatten;
		polys
		,
		{i,1,divFaces//Length}
	]//STIIMapleFindLROrders[#,vars,coeffs]&;
	commonLRorders=commonLRorders//SortBy[#,#[[2]]&]&;*)
	
	
	LRorders=Table[
		{polys,exps}=Iren[[i]]//STtoCoeffMonPols[#,vars]&//Last;
		polys=Table[If[Count[{exps[[i]]},eps,\[Infinity]]>0||exps[[i]]<0,polys[[i]],{}],{i,1,exps//Length}]//Flatten;
		(*(polys//dummySTIIMapleFindAllLROrdersEndingWith[#,vars,coeffs,Values[gauges][[i]]]&)//Print;*)
		LROrdersCompatibleWithGauges=(polys//STIIMapleFindAllLROrdersEndingWith[#,vars,coeffs,Complement[vars,#]&/@(Values[gauges][[i]])]&);
		Print[ToString[i]<>"/"<>ToString[(divFaces//Length)-1]<>"+1"];
		LROrdersCompatibleWithGauges//Length//Print;
		LROrdersCompatibleWithGauges=Join[LROrdersCompatibleWithGauges]//SortBy[#,#[[2]]&]&;
		If[Length[LROrdersCompatibleWithGauges]===0,Print["WARNING: no LR order compatible with gauge found for face: " <>ToString[i]]];
		LROrdersCompatibleWithGauges[[;;,1]];
		,
		{i,2,divFaces//Length}
	];
	
	(* Prepend LR orders for integrand *)
	Print["last one: ..."];
	LRorders=Join[{commonLRorders},LRorders];
	(*LRorders=Join[{(integrand//STtoCoeffMonPols[#,vars]&//Last//First//STIIMapleFindLROrders[{#},vars,coeffs]&)[[;;,1]]},LRorders];*)
	(*Return[LRorders];*)
	backup5={divFaces,Iren,LRorders,gauges//Values};
	
	(**)
	Join[{commonLRorders[[;;,1]]},Table[
	potentialLRorders=Intersection@@(LRorders[[Position[divFaces,Alternatives@@Cases[divFaces,v_/;SubsetQ[v,divFaces[[face]] ]],1]//Flatten]])[[;;,;;,;;(Length[vars]-Length[divFaces[[face]]])]];
	faceLRorders=potentialLRorders//Cases[#,v_/;Or@@(Table[Sort[v]===Sort[ga],{ga,Values[gauges][[face]]}])]&
	,{face,2,divFaces//Length}]]
	
];


*)


(* ::Subsection::Closed:: *)
(*Subtraction Formula  <- this is the one being used*)


(* Compute subtraction formula *) 
(* It is possible to pass tropical data and u-variables if these have already been found *)
Clear[STIISubtractionFormula];
Options[STIISubtractionFormula]={
	"Gauges"->"Automatic",
	"divUs"->{},
	"regulators"->{eps},
	"trDataGiven"->{}
};
STIISubtractionFormula[integrand_,vars_,coeffs_,OptionsPattern[]]:=Module[
	{
	divUs,
	regulators,
	trDataGiven,
	
	trData,regularizeFace,monspols,polsInt,divFacets,faces,divFaces,
	us,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,
	newIntegrand,result,A,B,trint,trValues,Wvects,dF,faceLROrders,storeLROrders},
	
	divUs=OptionValue["divUs"];
	regulators=OptionValue["regulators"];
	trDataGiven=OptionValue["trDataGiven"];
	
	Print["Computing tropical data"];
	trData=If[trDataGiven==={},integrand//STIntegrandTropicalData[#,vars,coeffs]&,trDataGiven];
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	If[Count[trValues/.Alternatives@@regulators->0,v_/;v> 0]>0,
		Print["Power divergence on "<>ToString[ Position[trValues/.Alternatives@@regulators->0,v_/;v> 0]//Flatten];]
	];
	
	
	monspols=STtoMonPols[integrand,vars];
	polsInt=Times@@(STtoMonPols[integrand,vars][[2,1]]);
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	(* If no divergences, return integrand in subtraction formula format.
	   Divide by Times@@vars to convert from dx/x to flat dx convention,
	   matching what regularizeFace does for non-trivial faces (line ~2507). *)
	If[divFacets==={},
		Print[" Integrand is already locally finite! " ];
		Return[{{{integrand/(Times@@vars)},vars}}];
	];
	
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	
	(*u-variables*)
	us=Table[{},{f,trData["rays"]}];
	If[Length[divUs]>0,
		Table[dF=divFacets[[j]];us[[dF]]=divUs[[j]],{j,1,divFacets//Length}]
	,
	Print["Computing u-variables "];
	(* STProduceUs[xvars_,rays_,divFaces_,asAnAssociationQ_:False,normalizeQ_:True,QuietQ_:True] *)
	us[[divFacets]]=If[divUs==={},
	(*STProduceAllUs[vars,trData["rays"],divFaces,True]//Factor//Print;*)
	STProduceUs[vars,trData["rays"],divFaces,False,True,True]//Factor
	,divUs//Factor];
	
	];
	If[Count[us,"NotFound!"]>0,
		Print[" The geometric property is not satisfied! "];
		Return[us];
	,
	Print["Proceeding with u-variables: " , us ]
	];
	
	
	omus=1-us//Factor;
	
	(* Recover the W-vectors from the u-vars *)
	Wvects= Table[
	A = (1-us[[dF]])//Factor//Numerator;
	B =( (1-us[[dF]])//Factor//Denominator)-(A//Factor)//Factor;
	dF-> CoefficientRules[B,vars][[1,1]]-CoefficientRules[A,vars][[1,1]]
	,{dF,divFacets}];
	Wvects=Association@@Wvects;
	(*J space*)
	(*Js=Table[
	If[fc==={1,19},
	Print["special"];
	fc-> {2,3,4,6,7,8,9},
	fc-> FirstCase[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
	]
	,{fc,divFaces[[2;;]]}
	];
	Js=Join[{divFaces[[1]]-> Range[vars//Length]},Js];
	Js=Association@@Js;*)
	Print["Choosing Face Gauges"];
	OptionValue["Gauges"]//Print;
	Switch[OptionValue["Gauges"],
		"Automatic",
		Print["Automatic Mode: Choosing any gauge"];
		Js=Table[
			fc-> FirstCase[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
			,
			{fc,divFaces[[2;;]]}
		];
		Js=Join[{divFaces[[1]]-> Range[vars//Length]},Js];
		Js=Association@@Js;
		backup2=Js;,
		
		"LR",
		Print["LR mode: Looking for LR face orders"];
		faceLROrders=STIIEnumerateFaceLROrders[integrand,vars,coeffs,us//DeleteCases[#,{}]&,{eps},trData];
		(*backup=faceLROrders;*)
		(*faceLROrders=backup;*)
		(*faceLROrders//Length//Print;*)
		storeLROrders=Association@@Table[divFaces[[j]]->faceLROrders[[j]],{j,1,divFaces//Length}];
		Js=Table[divFaces[[j]]->(Position[vars,Alternatives@@(faceLROrders[[j]][[1]])(* Arbitrary Choice *),1]//Flatten),{j,1,divFaces//Length}];
		Js=Association@@Js;
		backup3=Js;
	];
	Js//Print;
	(*Monomial*)
	extraMons=Association@@Table[
	fc->(1/(Times@@(vars^(-trValues[[fc]] . ( Wvects/@fc))))/.Alternatives@@vars[[Complement[Range[vars//Length],Js[fc]]]]->1)
	,{fc,divFaces}
	];
	
	(*Volume*)
	Print["Computing volumes"];
	volume=Table[
	fc->(1/Times@@(-trValues[[fc]]))(Abs[Det[Join[trData["rays"][[fc]],Table[KroneckerDelta[i,Js[fc][[j]]],{j,1,Length[Js[fc]]} ,{i,1,vars//Length}  ]]]])
	,{fc,divFaces[[2;;]]}
	];
	volume=Join[{divFaces[[1]]->1},volume];
	volume=Association@@volume;
	
	(*jacs = 1-u *)
	Print["Computing jacobians"];
	jacFace=Table[
	face->Table[Factor[(1-us[[j]])],{j,face}]
	,{face,divFaces}];
	jacFace = Association@@jacFace;
	
	(*We define a function which regularize a face*)
	regularizeFace[face_]:=Module[{compFacess,newIntegrandd},
		compFaces=Cases[divFaces,f_/;SubsetQ[f,face]];
		compFaces=Complement[#,face]&/@compFaces;
		(*We compute the face integrand to be regularized*)
		newIntegrand=If[face==={}, integrand,
		Times[
		STrestrictIntegrand[integrand,vars,trData["rays"][[face]]],
		extraMons[face]
		](*//STFactor//STSimplify*)
		];
		(*Print["computed new"];*)
		newIntegrand=newIntegrand;
		Monitor[
		result=volume[face] Table[
		If[subface==={},
		newIntegrand,
		Times[-a[subface],
		Times@@((jacFace[subface])^(1-trValues[[subface]]))(*why only 1+Trop, rather than n + Trop. Why set x[face]=0?*),
		STrestrictIntegrand[newIntegrand,vars,trData["rays"][[subface]]]
		]
		]
		,{subface,compFaces}];,
		{subface}];
		(*Print["computed result"];*)
		
		result=result(*//STFactor//STSimplify*);
		result=result/.a[subface_]:> (-1)^(Length[subface]+1);
		result=result/.Alternatives@@vars[[ Complement[Range[vars//Length],Js[face]] ]]->1;
		
		(*Print["computed gauge fix"];*)
		{result/(Times@@(vars[[Js[face]]])),Switch[OptionValue["Gauges"],"Automatic",vars[[Js[face]]],"LR",storeLROrders[face][[1]]]  }
	];
	Print["Computing subtraction formula"];
	Monitor[Table[regularizeFace[divFaces[[face]] ],{face,1,divFaces//Length}]//Return,{face,"/",divFaces//Length}]
];

(* Alias requested by Giulio: STTropicalSubtraction mirrors STSubtractionFormula. *)
STTropicalSubtraction = STSubtractionFormula;


(* ::Subsection:: *)
(*LR Subtraction Formula  Structure (no v_\rho's)*)


(* ::Subsubsection::Closed:: *)
(*Subtraction Formula Structure (with unspecified v's)*)


(* Compute subtraction formula *) 
(* It is possible to pass tropical data and u-variables if these have already been found *)
Clear[STIISubtractionFormulaStructure];
Options[STIISubtractionFormulaStructure]={
	"Gauges"->"Automatic",
	"regulators"->{eps},
	"trData"->{}
};
STIISubtractionFormulaStructure[integrand_,vars_,coeffs_,OptionsPattern[]]:=Module[
	{
	divUs,
	regulators,
	trDataGiven,
	
	trData,regularizeFace,monspols,polsInt,divFacets,faces,divFaces,
	us,omusPossibleValues,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,
	newIntegrand,result,A,B,trint,trValues,Wvects,dF,faceLROrders,storeLROrders},
	
	(*divUs=OptionValue["divUs"];*)
	regulators=OptionValue["regulators"];
	trDataGiven=OptionValue["trData"];
	
	Print["Computing tropical data"];
	trData=If[trDataGiven==={},integrand//STIntegrandTropicalData[#,vars,coeffs]&,trDataGiven];
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	If[Count[trValues/.Alternatives@@regulators->0,v_/;v> 0]>0,
		Print["Power divergence on "<>ToString[ Position[trValues/.Alternatives@@regulators->0,v_/;v> 0]//Flatten];]
	];
	
	
	monspols=STtoMonPols[integrand,vars];
	polsInt=Times@@(STtoMonPols[integrand,vars][[2,1]]);
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	(* If no divergences, return integrand in subtraction formula format.
	   Divide by Times@@vars to convert from dx/x to flat dx convention,
	   matching what regularizeFace does for non-trivial faces (line ~2507). *)
	If[divFacets==={},
		Print[" Integrand is already locally finite! " ];
		Return[{{{integrand/(Times@@vars)},vars}}];
	];
	
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	
	(*u-variables*)
	(*us=Table[v[f],{f,divFacets}];*)
	
	(*u-variables*)
	us=Table[{},{f,trData["rays"]}];
	If[Length[divUs]>0,
		Table[dF=divFacets[[j]];us[[dF]]=divUs[[j]],{j,1,divFacets//Length}]
	,
	Print["Computing u-variables "];
	(* STProduceUs[xvars_,rays_,divFaces_,asAnAssociationQ_:False,normalizeQ_:True,QuietQ_:True] *)
	us[[divFacets]]=STProduceAllUs[vars,trData["rays"],divFaces,True]//Factor;
	
	];
	If[Count[us,"NotFound!"]>0,
		Print[" The geometric property is not satisfied! "];
		Return[us];
	];
	
	omusPossibleValues=Association@@Table[f-> (1-us[[f]]//Factor),{f,divFacets}];
	omus=Association@@Table[f->v[f],{f,divFacets}];
	
	If[
		Count[us,"NotFound!"]>0,
		Print[" The geometric property is not satisfied! "];
		Return[us];
		,
		Print["Proceeding with v-functions: " , omus ]
	];
	
	(* Recover the W-vectors from the u-vars *)
	Wvects= Table[
	A = (1-us[[dF]])//Factor//Numerator;
	B =( (1-us[[dF]])//Factor//Denominator)-(A//Factor)//Factor;
	dF-> CoefficientRules[B,vars][[1,1]]-CoefficientRules[A,vars][[1,1]]
	,{dF,divFacets}];
	Wvects=Association@@Wvects;
	(*J space*)
	Print["Enumerating Face Gauges"];
	Js=Table[
			fc-> Cases[Subsets[vars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> Complement[Range[vars//Length],I]]
			,
			{fc,divFaces[[2;;]]}
		];
	Js=Join[{divFaces[[1]]-> {Range[vars//Length]}},Js];
	Js=Association@@Js;
	
	(* Warning: extra monomials, volumes and jacobians depend on a specific choice of u's. We simply omit these here*)
	(*Monomial*)
	(*extraMons=Association@@Table[
	fc->(1/(Times@@(vars^(-trValues[[fc]] . ( Wvects/@fc))))/.Alternatives@@vars[[Complement[Range[vars//Length],Js[fc]]]]->1)
	,{fc,divFaces}
	];*)
	extraMons=Association@@Table[
		fc->1
		,{fc,divFaces}
	];
	(*Volume*)
	
	volume=Table[
		fc->1
		,
		{fc,divFaces[[2;;]]}
	];
	volume=Join[{divFaces[[1]]->1},volume];
	volume=Association@@volume;
	
	jacFace=Table[
	face->Table[omus[j],{j,face}]
	,{face,divFaces}];
	jacFace = Association@@jacFace;
	(*
	jacFace=Table[
	face->Table[Factor[(1-us[[j]])],{j,face}]
	,{face,divFaces}];
	jacFace = Association@@jacFace;*)
	
	(*We define a function which regularize a face*)
	regularizeFace[face_]:=Module[{compFacess,newIntegrandd},
		compFaces=Cases[divFaces,f_/;SubsetQ[f,face]];
		compFaces=Complement[#,face]&/@compFaces;
		(*We compute the face integrand to be regularized*)
		newIntegrand=If[face==={}, 
			integrand
			,
			Times[
				STrestrictIntegrand[integrand,vars,trData["rays"][[face]]],
				extraMons[face]
				]
			];
		(*Print["computed new"];*)
		newIntegrand=newIntegrand;
		Monitor[
			result=volume[face] Table[
				If[subface==={},
					newIntegrand,
					Times[-a[subface],
						Times@@((jacFace[subface])^(1-trValues[[subface]]))(*why only 1+Trop, rather than n + Trop. Why set x[face]=0?*),
						STrestrictIntegrand[newIntegrand,vars,trData["rays"][[subface]]]
					]
				]
				,
				{subface,compFaces}
			];
			,
			{subface}
		];
		(*Print["computed result"];*)
		
		result=result(*//STFactor//STSimplify*);
		result=result/.a[subface_]:> (-1)^(Length[subface]+1);
		(* We don't chose a gauge at this stage *)
		(*result=result/.Alternatives@@vars[[ Complement[Range[vars//Length],Js[face]] ]]->1;*)
		
		(*Print["computed gauge fix"];*)
		{result,vars[[#]]&/@(Js[face]) }
	];
	Print["Computing subtraction formula"];
	Monitor[
	{divFaces,Table[regularizeFace[divFaces[[face]] ],{face,1,divFaces//Length}],omusPossibleValues}//Return
	
	,{face,"/",divFaces//Length}]
];


(* ::Subsubsection::Closed:: *)
(*Scan Over v's and Gauges*)


Clear[STIIFindLRvs];
STIIFindLRvs[subtractionStructure_,xvars_,coeffs_]:=Module[
{
divFaces,
divFacets,
choice,
face,
varsForThisProduct,
choicesOfVs,
checkQ,
Iren,listsOfPolys,
faceLRorders,
LRorders=<||>,
LRchoice,
polys,exponents
},

divFaces=subtractionStructure[[1]];
divFacets=subtractionStructure[[-1]]//Keys;
choicesOfVs=Tuples[subtractionStructure[[-1]]//Values];
(*choicesOfVs=Join[{1-{x3/(x2+x3),x7/(x7+x8),(x10 x7)/(x10 x7+x6 x9),x10/(1+x10),x6/(x6+x7),x5/(1+x5),1/(1+x10),x1/(x1+x2)}//Factor},choicesOfVs];*)
choicesOfVs=Table[Association@@(Table[divFacets[[i]]-> choice[[i]],{i,1,divFacets//Length}]),{choice,choicesOfVs}];
choicesOfVs[[;;2]]//Print;
(* loop over choices of vs *)
Do[
	choice=choicesOfVs[[k]];
	Print[k,"/",choicesOfVs//Length," : ", choice];
	
	checkQ=True;

	(*(* Faster pre-check*)
	Do[
		face=divFaces[[ii]];
		If[Length[face]===1,Continue[] (* a single v is always LR*)];
		(*check if vs of this face are LR *)
		varsForThisProduct=(Variables/@Denominator/@choice/@face)//Flatten//DeleteDuplicates;
		(*wrapSTIIMapleFindLROrders[{(Denominator/@choice/@face)},varsForThisProduct,coeffs]//Print;*)
		checkQ= Length[STIIMapleFindLROrders[{(Denominator/@choice/@face)},varsForThisProduct,coeffs]]>0;
		If[checkQ===False,If[checkQ===False,Print["pre-check failed on face ",face," for choice ", choice];Break[]];,Print["pre-check passed on this face"]]
		,
		{ii,2,divFaces//Length}
	];
	If[checkQ===False,Break[]];*)
	
	(* Loop over faces *)
	Do[
	Print["face: ",i];
		(* loop over choices of gauges *)
		Do[
			(*DeleteCases[xvars,Alternatives@@gauge]//Print["gauge",#]&;*)
			Iren=subtractionStructure[[2,i,1]]/.v[a_]:>choice[a]/.Table[v->1,{v,DeleteCases[xvars,Alternatives@@gauge]}];
			listsOfPolys=Table[
				{polys,exponents}=(STtoCoeffMonPols[ct,gauge])//Last;
				polys[[Position[exponents,v_/;Count[{v},eps,\[Infinity]]>0 || v<0,{1}]//Flatten]]
			,{ct,Iren}];
			(*listsOfPolys//Print;*)
			faceLRorders=STIIMapleFindLROrders[listsOfPolys,gauge,coeffs];
			(*{Table[v->1,{v,DeleteCases[xvars,Alternatives@@gauge]}],Iren,wrapSTIIMapleFindLROrders[listsOfPolys,gauge,coeffs],gauge,DeleteCases[xvars,Alternatives@@gauge]}//Print;*)
			checkQ= Length[faceLRorders]>0;
			(* We stop as soon as we find a LR gauge*)
			If[checkQ===True,LRorders[divFaces[[i]]]=faceLRorders;Break[]]
			,
			{gauge,subtractionStructure[[2,i,2]]}
		];
		(* We stop as soon as we find a face with no LR gauges *)
		If[checkQ===False, Print["this face has no LR gauges!"]; Break[]]
		,
		{i,1,Min[\[Infinity],divFaces//Length]}
	];
	If[checkQ===True,Break[]]
	,
	{k,1,choicesOfVs//Length}
];
If[checkQ===False,
	Print["WARNING: no choice of v-functions is LR!"];
	Return[{{},{}}],
	Print["Success!"]; 
	Return[{choice,LRorders}]
];

]


(* ::Subsubsection:: *)
(*Subtraction Formula With Given v's and Gauges*)


(* Compute subtraction formula *) 
(* It is possible to pass tropical data and u-variables if these have already been found *)
Clear[STIISubtractionFormulaFromScan];
Options[STIISubtractionFormulaFromScan]={
	"Gauges"->"Automatic",
	"regulators"->{eps},
	"trData"->{}
};
STIISubtractionFormulaFromScan[integrand_,vars_,coeffs_,scan_,OptionsPattern[]]:=Module[
	{trData,regularizeFace,monspols,polsInt,divFacets,faces,divFaces,
	us,omus,vvWvects,Js,extraMons,jacFace,volume,compFaces,
	newIntegrand,result,A,B,trint,trValues,Wvects},
	
	(*divUs=OptionValue["divUs"];*)
	regulators=OptionValue["regulators"];
	trDataGiven=OptionValue["trData"];
	
	
	Print["Computing tropical data"];
	trData=If[trDataGiven==={},integrand//STIntegrandTropicalData[#,vars,coeffs]&,trDataGiven];
	
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;
	If[Count[trValues/.Alternatives@@regulators->0,v_/;v> 0]>0,
		Print["Power divergence on "<>ToString[ Position[trValues/.Alternatives@@regulators->0,v_/;v> 0]//Flatten];]
	];
	
	monspols=STtoMonPols[integrand,vars];
	polsInt=Times@@(STtoMonPols[integrand,vars][[2,1]]);
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	(* If no divergences, return integrand in subtraction formula format.
	   Divide by Times@@vars to convert from dx/x to flat dx convention,
	   matching what regularizeFace does for non-trivial faces (line ~2507). *)
	If[divFacets==={},
		Print[" Integrand is already locally finite! " ];
		Return[{{{integrand/(Times@@vars)},vars}}];
	];
	
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	
	(*u-variables*)
	us=scan[[1]];
	(*Return[];*)
	(* Note: this syntax works also on Associations *)
	omus=1-us//Factor;
	
	(* Recover the W-vectors from the u-vars *)
	Wvects= Table[
	A = (1-us[dF])//Factor//Numerator;
	B =( (1-us[dF])//Factor//Denominator)-(A//Factor)//Factor;
	dF-> CoefficientRules[B,vars][[1,1]]-CoefficientRules[A,vars][[1,1]]
	,{dF,divFacets}];
	Wvects=Association@@Wvects;
	(*J space*)
	Js=First/@Last/@scan[[2]];
	(* Scan lists the actual variables, we convert it into positions in the list vars*)
	Js=(Flatten[Position[vars, #] & /@ #]&)/@Js;
	(*Monomial*)
	extraMons=Association@@Table[
	fc->(1/(Times@@(vars^(-trValues[[fc]] . ( Wvects/@fc))))/.Alternatives@@vars[[Complement[Range[vars//Length],Js[fc]]]]->1)
	,{fc,divFaces}
	];
	
	(*Volume*)
	Print["Computing volumes"];
	volume=Table[
	fc->(1/Times@@(-trValues[[fc]]))(Abs[Det[Join[trData["rays"][[fc]],Table[KroneckerDelta[i,Js[fc][[j]]],{j,1,Length[Js[fc]]} ,{i,1,vars//Length}  ]]]])
	,{fc,divFaces[[2;;]]}
	];
	volume=Join[{divFaces[[1]]->1},volume];
	volume=Association@@volume;
	
	(*jacs = 1-u *)
	Print["Computing jacobians"];
	jacFace=Table[
	face->Table[Factor[(1-us[j])],{j,face}]
	,{face,divFaces}];
	jacFace = Association@@jacFace;
	
	(*We define a function which regularize a face*)
	regularizeFace[face_]:=Module[{compFacess,newIntegrandd},
	compFaces=Cases[divFaces,f_/;SubsetQ[f,face]];
	compFaces=Complement[#,face]&/@compFaces;
	(*We compute the face integrand to be regularized*)
	newIntegrand=If[face==={}, integrand,
	Times[
	STrestrictIntegrand[integrand,vars,trData["rays"][[face]]],
	extraMons[face]
	](*//STFactor//STSimplify*)
	];
	(*Print["computed new"];*)
	newIntegrand=newIntegrand;
	Monitor[
	result=volume[face] Table[
	If[subface==={},
	newIntegrand,
	Times[-a[subface],
	Times@@((jacFace[subface])^(1-trValues[[subface]]))(*why only 1+Trop, rather than n + Trop. Why set x[face]=0?*),
	STrestrictIntegrand[newIntegrand,vars,trData["rays"][[subface]]]
	]
	]
	,{subface,compFaces}];,
	{subface}];
	(*Print["computed result"];*)
	
	result=result(*//STFactor//STSimplify*);
	result=result/.a[subface_]:> (-1)^(Length[subface]+1);
	result=result/.Alternatives@@vars[[ Complement[Range[vars//Length],Js[face]] ]]->1;
	(*Print["computed gauge fix"];*)
	
	{result/(Times@@(vars[[Js[face]]])),vars[[Js[face]]]}
	];
	Print["Computing subtraction formula"];
	Monitor[Table[regularizeFace[divFaces[[face]] ],{face,1,divFaces//Length}]//Return,{face,"/",divFaces//Length}]
];

(* Alias requested by Giulio: STTropicalSubtraction mirrors STSubtractionFormula. *)
STTropicalSubtraction = STSubtractionFormula;


(* ::Subsection::Closed:: *)
(*Find NP continuation*)


(* Find the simplest NP continuation *)
Clear[STIIfindFirstNPContinuation];
STIIfindFirstNPContinuation[vars_,trData_,divFacets_,startFrom_:0,forcedUs_:{}]:=Module[
	{
	foundQ,hypus,howManyRaysContinue,continuedRays,result,counter,subs
	},
	
	foundQ=False;
	(* If imposed to use some u's, no NP supposed to be done *)
	foundQ=Not[forcedUs==={}];
	If[foundQ,Return[{{},forcedUs}]];
	
	howManyRaysContinue=startFrom;

	While[Not[foundQ],
		Print["Assuming "<>ToString[howManyRaysContinue]<>" divergent rays have been deleted by NP..."];
		subs=Subsets[divFacets,{howManyRaysContinue}];
		counter=1;
		While[Not[foundQ] &&(counter<= Length[subs]),
			(* Assuming these rays are no longer there... *)
			continuedRays=subs[[counter]];
			(* ... we would find these hypotetical u's: *)
			hypus= STProduceUs[vars,trData["rays"],STGetFaces[trData,Complement[divFacets,continuedRays]],True,True];
			(* Check if we found satisfactory u's *) 
			(*forbidden={1}|{8}|{12}|{13}|{25}|{30};*)
			forbidden={2}|{7};
			If[
				Count[hypus,"NotFound"|(Alternatives@@vars)^(n_/;n>1),\[Infinity]]>0 || Count[Denominator/@hypus,den_/;Count[den,Alternatives@@vars,\[Infinity]]>(multiDegreeNP)]>0||MatchQ[subs[[counter]],forbidden],
				foundQ=False,
				Print["... found!"];
				foundQ=True;
				STProduceAllUs[vars,trData["rays"],STGetFaces[trData,Complement[divFacets,continuedRays]],True]//Print;
				Break[];
			];
			counter++;
		];
		howManyRaysContinue++;
		If[howManyRaysContinue>Length[divFacets],
		foundQ=False; 
		Break[]
		];
	];
	If[foundQ,
		Return[{subs[[counter]],hypus}],
		"STfindFirstNPContinuation::Something went wrong, we have not found a continuation!"//Print;
		Return[{NotFound,NotFound}]
	];
];


(* ::Subsection::Closed:: *)
(*Main Method: Subtraction+NP*)


(* Computes tropical subtraction formula + NP continuation *) 
(* It looks for the simplest NP continuation reducing to integrands with the geometric property and applies the subtraction formula to it *)
(* It is possible to pass tropical data and u-variables if these have already been found *)
Clear[STIIExpandIntegral];
Options[STIIExpandIntegral]={
	"Gauges" -> "Automatic",
	"regulators"->{eps},
	"trDataGiven"-> {},
	"forceUs"-> {}
};
STIIExpandIntegral[integrand_,vars_,coeffs_,opts:OptionsPattern[]]:=Module[
	{
	regulators,trDataGiven,forceUs,
	(* prefactor *)
	pref,mons,pols, intnopref,
	(* Tropical Analysis *)
	trData,trint,trValues,divFacets,faces,divFaces,
	wherePowerDiv,ordersOfDivergences,raysToContinue,
	(* Subtraction Formula*)
	us,omus,vvWvects,Wvects,Js,extraMons,jacFace,volume,
	regularizeFace,compFaces,newIntegrand,
	result,A,B,
	(* NP reduction *)
	hypus,
	NPcontinuationSub, NPcontinuationPow, NPcontinuation, ordersPref, subtractableInts,
	expansions (* results *)
	}
	,
	regulators=OptionValue["regulators"];
	trDataGiven=OptionValue["trDataGiven"];
	forceUs=OptionValue["forceUs"];
	
	

	(* Zero-variable short-circuit: an empty integration domain is trivially
	   locally finite, with the whole integrand playing the role of prefactor.
	   Return it in the same shape as the "Locally finite integrand" branch
	   below so downstream consumers can process it uniformly. *)
	If[vars === {},
		Return[{{integrand, {{{1}, {}}}}}]
	];

	(* Extract prefactor, monomials and polynomials *)
	{pref,mons,pols}=STtoCoeffMonPols[integrand,vars];
	intnopref=(Times@@(vars^mons))(Times@@(pols[[1]]^pols[[2]]));
	Print[pref, " ", intnopref ];
		
	(* Computing Tropical Data *)
	Print["Loading or Computing tropical data"];
	If[trDataGiven==={},
		trData=integrand//STIntegrandTropicalData[#,vars,coeffs]&,
		trData=trDataGiven;
	];
	trData//Print;
	
	trint=integrand//STTropicalizeIntegrand[#,vars,coeffs]&;
	trValues=Table[STEvalRay[trint,ray,vars],{ray,trData["rays"]}]//Factor;

	(* Check if trop = 0 on some ray*)
	If[MemberQ[trValues,0],
		Print[" The Euler integral is not well defined! Aborting ... "];
		(* PROBE: Message survives runQuiet's Print->Null rebinding; if you see
		   STExpandIntegral::notwelldef next run, the guard is confirmed as the
		   cause of the upstream $Aborted. Then swap Abort[] -> Return[$Aborted]
		   for the real fix. *)
		Message[STExpandIntegral::notwelldef];
		Abort[];
	];
	
	
	divFacets=Position[trValues/.Alternatives@@regulators->0,v_/;v>= 0]//Flatten;
	faces=trData//STGetFaces[#,divFacets]&;
	divFaces=Cases[faces,v_/;SubsetQ[divFacets,v]];
	ordersOfDivergences=Association@@Table[facet-> (trValues[[facet]]/.eps->0),{facet,divFacets}];
	
	(* If already locally finite nothing to do *)
	If[divFacets==={},
		Print["Locally finite integrand, no need for subtractions!"];
		Return[{{pref,{{{intnopref/(Times@@vars)}, vars}}}}]
	];		
	
	If[MemberQ[Values[ordersOfDivergences],v_/;v>0],
		wherePowerDiv=divFacets[[ Position[Values[ordersOfDivergences],v_/;v>0]//Flatten ]];
		Print["Power divergences on rays "<>ToString[wherePowerDiv]];
		Print["We will later apply Nillson-Passare continuation along these rays!"],
		wherePowerDiv={}
	];
	(* Determining minimal Nillson-Passare continuation *) 
	Print["Determining minimal Nillson-Passare continuation!"];
	{raysToContinue,hypus}=STIIfindFirstNPContinuation[vars,trData,divFacets,0,forceUs];
	Print["To apply formula, Nillson-Passare continue along rays: "<>ToString[raysToContinue]];
	
	
	(* Computing minimal Nillson-Passare continuation *) 
	Print["Computing minimal Nillson-Passare continuation!"];
	(* How many times each ray? *)
	NPcontinuationSub=Join@@Table[
		(* These must be made finite *)
		Table[trData["rays"][[ray]],ordersOfDivergences[ray]+1]
	,
	{ray,raysToContinue}
	];
	NPcontinuationPow=Join@@Table[
			(* These must be made log-divergent *)
			Table[trData["rays"][[ray]],ordersOfDivergences[ray]]
		,
		{ray,Complement[wherePowerDiv,raysToContinue]}
	];
	NPcontinuation=Join[NPcontinuationSub,NPcontinuationPow];
	
	(* Compute continuations *)
	If[NPcontinuation==={},
		subtractableInts={{pref,intnopref}};
		ordersPref={0};
		,
		subtractableInts=STContinueRays[{{pref,intnopref}},vars,NPcontinuation];
		ordersPref=subtractableInts[[;;,1]]dummy^eps//Series[#,{eps,0,0}]&;
		ordersPref=ordersPref[[;;,4]];
	];
	Print["----------"];
	Print["The reduction produced "<>(subtractableInts//Length//ToString)<>" integrands with the geometrical property!"];
	Print["They have prefactors of orders: "<>(eps^ordersPref//ToString[#,InputForm]&)];
	Print["----------"];
	
	(* Computing Subtraction Formulae *)
	Print["Computing Subtractions for each integrand with the geometrical property"];
	Print["----------"];
	expansions=Table[
	{int[[1]],STIISubtractionFormula[int[[2]], vars, coeffs,opts]}
	,{int,subtractableInts[[;;1]]}
	];
	Print["----------"];
	expansions
];


(* ::Subsection::Closed:: *)
(*Ancillary (Integrate Expansion with Mathematica)*)


(*(* Integrate a subtraction directly in Mathematica *)
(*STIntegrateSubtractionNP[expansion_,order_,integrationFunction_:Integrate, options_: Sequence[] ]:=
Table[
	Assuming[
		STPositiveVariables,
		integrationFunction[
			exp[[1]]//SeriesCoefficient[#,{eps,0,order}]&//Total//Factor//Evaluate
			,
			Evaluate[Sequence@@({#,0,\[Infinity]}&/@exp[[2]])], options
		]
	]
	,{exp,expansion}]//Total
];*)

(* Integrate a single subtraction directly in Mathematica *)
STIntegrateSingleSubtraction[expansion_, order_, integrationFunction_:Integrate, options_:{} ]:=Module[
{},
	Table[
		integrationFunction[ 
			integrand[[1]]//SeriesCoefficient[#,{eps,0,order}]&//Total , 
			integrand[[2]] , Sequence@@options
		]
		,
		{integrand, expansion}
	]//Total	
]*)


(*(* Integrate subtraction(s) directly in Mathematica *)
STIntegrateSubtractionNP[expansions_,order_,integrationFunction_:Integrate, options_: {} ]:=Module[
	{
	minOrder,maxOrder,
	results
	}
	,
	results=Table[
		minOrder=(expansions[[cnt,2,-1,2]]//Length)-(expansions[[cnt,2,1,2]]//Length);
		maxOrder=order + Max[0,-Series[expansions[[cnt,1]]dummy^eps,{eps,0,0}][[4]]];
		Print["Integrating contribution : "<>ToString[cnt]<>" from eps-order "<>ToString[minOrder] <>" to "<>ToString[maxOrder]];
		{
			expansions[[cnt,1]],
			Table[
				Print["integrating order: "<>ToString[ord]];
				{eps^ord , STIntegrateSingleSubtraction[expansions[[cnt,2]], ord, integrationFunction, options ]}
			,
			{ord, minOrder, maxOrder}
			]
		}
		,
		{cnt,1,expansions//Length}
	];
	(* Combine Results *)
	STCombineResultsMaple[results,order]
];*)


(*STIntegrateOrders[integrand_, vars_, maxTime_:60 , options_:{}]:=Module[
	{foundQ,currentMaxtime,timeIncrement,outTimeConst}
	,
	currentMaxtime=maxTime;
	timeIncrement=10;
	Do[
		outTimeConst=TimeConstrained[
				Assuming[
					STPositiveVariables
					,
					Integrate[integrand , Evaluate@(Sequence@@({#,0,\[Infinity]}&/@ perm)), Evaluate[Sequence@@options], GenerateConditions->False ]
				]
				,
				currentMaxtime,"IntegrationFailed"
			];
		foundQ=Not["IntegrationFailed"===outTimeConst];
		If[foundQ,outTimeConst//Print; Return[outTimeConst],"Failed"]
		,
		{perm,Permutations[vars//RandomSample]}
	]
];
*)


(*Clear[STFastIntegration];
STFastIntegration[integrand_,xvars_, timeLimit_:{1}, options_:{}]:=
Module[
{
nextVar,remainingVars,result,
tryIntegration,tryAllPermutations,success=False,value,
partialIntegrand,
MaxTimeLimit,timeLimitStep,currentTimeLimit,newPartialIntegrand,pPrint
},
	timeLimitStep=5;
	MaxTimeLimit=30;
	currentTimeLimit=timeLimit;
	
	Print["Current time limit is: "<>stFormatTime[currentTimeLimit]];
	
	tryIntegration[int_,var_]:=TimeConstrained[Assuming[STPositiveVariables,Integrate[int, {var,0,\[Infinity]},GenerateConditions->False]],currentTimeLimit,$Failed];
	
	remainingVars=xvars;
	partialIntegrand=integrand;
	
	While[(currentTimeLimit<=MaxTimeLimit)&&Not[success],
		Print["Trying time limit: "<>stFormatTime[currentTimeLimit]];
		Print["Remaining variables: "<>ToString[remainingVars]];
		Print["Current Integrand: "<>ToString[partialIntegrand,InputForm]];
		
		remainingVars=remainingVars//RandomSample;
		If[Length[remainingVars]===1,currentTimeLimit=MaxTimeLimit];
		Do[
			Print["trying ", nextVar];
			newPartialIntegrand=tryIntegration[partialIntegrand,nextVar];
			If[newPartialIntegrand=!=$Failed,
				Print["worked!",newPartialIntegrand];
				currentTimeLimit=timeLimit;
				remainingVars=DeleteCases[remainingVars, nextVar];
				partialIntegrand=newPartialIntegrand;
				Break[];
			]
			,
			{nextVar,remainingVars} 
		];
		success=Length[remainingVars]===0;
		If[newPartialIntegrand===$Failed,currentTimeLimit+=timeLimitStep;];
	];
	If[success,partialIntegrand,$Failed]
]*)


(*Clear[STFastIntegration2];
STFastIntegration2[integrand_,xvars_, timeLimit_:{1}, options_:{}]:=
Module[
{
nextVar,remainingVars,result,
tryIntegration,tryAllPermutations,success=False,value,
partialIntegrand,
MaxTimeLimit,timeLimitStep,currentTimeLimit,newPartialIntegrand,pPrint,
sortingFunction
},
	timeLimitStep=30;
	MaxTimeLimit=600;
	currentTimeLimit=timeLimit;
	
	Print["Integrating ", integrand, " over ", xvars];
	remainingVars=xvars;
	partialIntegrand=integrand;
	sortingFunction[v_]:=ByteCount[v]+ \[Infinity] Count[v,(HypergeometricPFQ|Hypergeometric2F1|Hypergeometric1F1)[___],\[Infinity]];
	
	While[(currentTimeLimit<MaxTimeLimit)&&(Length[remainingVars]>0),
	
		If[Length[remainingVars]===1, currentTimeLimit=MaxTimeLimit];
		Print["Current time limit is: "<>stFormatTime[currentTimeLimit]];
		Print["Remaining variables: "<>ToString[remainingVars]];
		Print["Current Integrand: "<>ToString[partialIntegrand,InputForm]];
		Print[" "];
		
		newPartialIntegrand=Table[
			{v,TimeConstrained[Assuming[STPositiveVariables,Integrate[ partialIntegrand , {v,0,\[Infinity]} , GenerateConditions->False ]],currentTimeLimit,"Failed"]}
			,
			{v,remainingVars}
		];
		newPartialIntegrand=DeleteCases[newPartialIntegrand,{_,"Failed"}];S
		If[newPartialIntegrand==={}, 
			currentTimeLimit=currentTimeLimit+timeLimitStep;
			Print["No variable was integrated, increasing allowed time..."];
			,
			newPartialIntegrand=SortBy[newPartialIntegrand,(#[[2]]//sortingFunction)&];
			partialIntegrand=newPartialIntegrand[[1,2]];
			remainingVars=DeleteCases[remainingVars, newPartialIntegrand[[1,1]] ]
		];
	];
	success=Length[remainingVars]===0;
	(*Print["done!", partialIntegrand, remainingVars];*)
	
	If[success, Print["Success!"]; partialIntegrand, Print["Failure!"]; $Failed ]
]

*)


(* ::Section::Closed:: *)
(*Maple Interface*)


(* ::Subsection:: *)
(*Ask Maple  To...*)


(* ::Subsubsection::Closed:: *)
(*Integrate*)


(* Ask Maple to integrate a single integrand with variables from 0 to infinity *)
STIIMapleHyperInt[integrand_,vars_,fibrationBasisQ_:False, rangeMin_:"0",rangeMax_:"infinity"]:=Module[{toIntegrate,domain,prefix,suffix,integrationString,saveString,fileMapleCommand,fileMapleOut,flatIntegrand,flatVars,unflatten,rawResult},

(* Flatten any subscripted atoms (MM[1], mm[2], ...) to flat names so
   Maple's FromMma translator doesn't misread them as function calls.
   Apply jointly to integrand + vars so the substitution is consistent;
   inverse rules are applied to the returned expression. *)
{flatIntegrand, unflatten} = SubTropica`stFlattenIndexedSymbols[{integrand, vars}];
flatVars = flatIntegrand[[2]];
flatIntegrand = flatIntegrand[[1]];

(* Integration String *)

domain=(StringJoin@@Table[ToString[v]<>"="<>rangeMin<>".."<>rangeMax<>",",{v,flatVars}])<>"]";
domain="["<>StringReplace[domain,",]"-> "]"];
toIntegrate=flatIntegrand//ToString[#,InputForm]&//SubTropica`stStripMmaContexts;
prefix="with(MmaTranslator);\n read \""<>SubTropica`$SThyperIntPath<>"\";\n _hyper_abort_on_divergence :=false;\n result:=hyperInt(FromMma(`";
suffix="`),"<>domain<>");";
integrationString=prefix<>toIntegrate<>suffix;
If[fibrationBasisQ,
integrationString=integrationString<>"\n finalRes:=fibrationBasis(result):",
integrationString=integrationString<>"\n finalRes:=result:"
];


(* Result String*)
fileMapleOut=CreateFile[];
saveString="fileRes:=\""<>fileMapleOut<>"\":\n";
saveString=saveString<>"if FileTools[Exists](fileRes) then FileTools[Remove](fileRes): end if:\nFileTools[Text][Open](fileRes):\nFileTools[Text][WriteString](fileRes,convert(finalRes,string)):\nFileTools[Text][Close](fileRes):";

fileMapleCommand=CreateFile[];
WriteString[fileMapleCommand,integrationString<>"\n"<>saveString];
Close[fileMapleCommand];

(* Launch Maple*)
(*fileMapleCommand//Print;*)
Run[SubTropica`$MapleCommand<>" -q "<>fileMapleCommand];
rawResult = If[fibrationBasisQ,
ReadString[fileMapleOut]//SubTropica`STEvaluateII,
(* Chain guard: if STformatHyperIntMapleOut returned $Failed (Maple output
   contained `Root(...)`, not round-trippable), propagate $Failed instead
   of letting EvaluatePeriods quietly wrap a broken input. *)
Module[{parsed=ReadString[fileMapleOut]//SubTropica`STformatHyperIntMapleOut},
    If[parsed===$Failed,$Failed,HyperIntica`EvaluatePeriods[parsed]]]
];
(* Unflatten: restore MM__1 -> MM[1], etc. *)
If[rawResult === $Failed, $Failed, rawResult /. unflatten]

];


(* ::Subsubsection::Closed:: *)
(*Find LR Orders*)


STIIMapleFindLROrders[listsOfpolys_,xvars_,coeffs_]:=Module[{FirstScoredOrders,result},
	FirstScoredOrders=STIIMapleFindAllLROrders[listsOfpolys//First,xvars,coeffs];
	STIIMapleScoreOrders[FirstScoredOrders,Drop[listsOfpolys,{1}],xvars,coeffs]//SortBy[#,-Last[#]&]&
]


STIIMapleScoreOrders[ScoredOrders_,listsOfpolys_,xvars_,coeffs_]:=Module[{newScores},
	(*Print[ScoredOrders];*)
	(* New Scores *)
	newScores=STIIMapleScoreOrdersForPolys[ScoredOrders[[;;,1]],listsOfpolys//First,xvars,coeffs];
	(* Add Old Scores *)
	newScores=Table[nsc[[1]]-> nsc[[2]]+(nsc[[1]]/.ScoredOrders),{nsc,newScores}];
	(* Iterate *)
	STIIMapleScoreOrders[newScores,Drop[listsOfpolys,1],xvars,coeffs]
]


STIIMapleScoreOrders[ScoredOrders_,{},xvars_,coeffs_]:=ScoredOrders


(* Ask Maple to keep only LR orders and score them  *)
Clear[STIIMapleScoreOrdersForPolys];
STIIMapleScoreOrdersForPolys[orders_,polys_,xvars_,coeffs_]:=Module[{
polysAndPairs,
prefix,saveString,fileMapleCommand,
fileMapleOut,
rawResult
},

	
	polysAndPairs="["<>({polys,Subsets[polys,{2}]}//ToString[#,InputForm]&//stStripMmaContexts//StringTake[#,2;;-2]&)<>"]:";
	

	(* Main Maple Script *)
	prefix="with(MmaTranslator);\nread \""<>$SThyperIntPath<>"\";\n";
	prefix=prefix<>"vars:=["<>(ToString[xvars]//StringTake[#,2;;-2]&)<>"]:\n";

	prefix=prefix<>("unassign(`L`):L[{}]:="<>polysAndPairs<>"\n");
	prefix=prefix<>"cgReduction(L,"<>stStripMmaContexts[ToString[coeffs]]<>"):\n";
	
	prefix=prefix<>STIIMapleModCheckIntegrationOrder<>"\n"<>STIIMapleModFindAllOrders<>"\n"<>STIIMapleModFindOrdersWithScore<>"\n";
	prefix=prefix<>STIIMapleModScoreOrders<>"\n";
	prefix=prefix<>"orders:="<>(orders//ToString//StringReplace[#,{"{"->"[","}"->"]"}]&)<>":\n";
	
	
	
	prefix=prefix<>"finalRes:=ModScoreOrders(orders,L):\n";
	(*prefix//Print;*)
	
	(* Result String*)
	fileMapleOut=CreateFile[];
	saveString="fileRes:=\""<>fileMapleOut<>"\":\n";
	saveString=saveString<>"if FileTools[Exists](fileRes) then FileTools[Remove](fileRes): end if:\n FileTools[Text][Open](fileRes):\n";
	saveString=saveString<>"FileTools[Text][WriteString](fileRes,convert(finalRes,string)):\n";
	saveString=saveString<>"FileTools[Text][Close](fileRes):";
	
	(* Save Maple Command *)
	fileMapleCommand=CreateFile[];
	(*prefix<>"\n"<>saveString//Print;*)
	WriteString[fileMapleCommand,prefix<>"\n"<>saveString];
	Close[fileMapleCommand];
	
	Run[$MapleCommand<>" -q "<>fileMapleCommand];
	(*ReadString[fileMapleOut]//Print;*)
	ReadString[fileMapleOut]//StringReplace[#,{"\""->"",  "["->"{",  "]"-> "}","="-> "->" }]&//ToExpression

]


(* Ask Maple to compute all LR orders for a list of polynomials *)
Clear[STIIMapleFindAllLROrders];
STIIMapleFindAllLROrders[polys_,xvars_,coeffs_]:=Module[{
polysAndPairs,
prefix,saveString,fileMapleCommand,
fileMapleOut,
rawResult
},

	polysAndPairs="["<>({polys,Subsets[polys,{2}]}//ToString[#,InputForm]&//stStripMmaContexts//StringTake[#,2;;-2]&)<>"]:";

	(* Main Maple Script *)
	prefix="with(MmaTranslator);\nread \""<>$SThyperIntPath<>"\";\n";
	prefix=prefix<>"vars:=["<>(ToString[xvars]//StringTake[#,2;;-2]&)<>"]:\n";

	prefix=prefix<>("unassign(`L`):L[{}]:="<>polysAndPairs<>"\n");
	prefix=prefix<>"cgReduction(L,"<>stStripMmaContexts[ToString[coeffs]]<>"):\n";
	prefix=prefix<>STIIMapleModCheckIntegrationOrder<>"\n"<>STIIMapleModFindAllOrders<>"\n"<>STIIMapleModFindOrdersWithScore<>"\n";
	
	prefix=prefix<>"orders:=ModFindAllOrders(L,vars):\n";
	prefix=prefix<>"finalRes:=orders:\n";
	(*prefix//Print;*)
	
	(* Result String*)
	fileMapleOut=CreateFile[];
	saveString="fileRes:=\""<>fileMapleOut<>"\":\n";
	saveString=saveString<>"if FileTools[Exists](fileRes) then FileTools[Remove](fileRes): end if:\n FileTools[Text][Open](fileRes):\n";
	saveString=saveString<>"FileTools[Text][WriteString](fileRes,convert(finalRes,string)):\n";
	saveString=saveString<>"FileTools[Text][Close](fileRes):";
	
	(* Save Maple Command *)
	fileMapleCommand=CreateFile[];
	WriteString[fileMapleCommand,prefix<>"\n"<>saveString];
	Close[fileMapleCommand];
	
	Run[$MapleCommand<>" -q "<>fileMapleCommand];
	(*ReadString[fileMapleOut]//Print;*)
	ReadString[fileMapleOut]//StringReplace[#,{"\""->"",  "["->"{",  "]"-> "}","="-> "->" }]&//ToExpression

]


(* ::Subsubsection::Closed:: *)
(*Find LR Orders: Linear Crawl*)


STIIMapleFindLROrders[listsOfpolys_,xvars_,coeffs_]:=Module[{FirstScoredOrders},
	FirstScoredOrders=STIIMapleFindAllLROrders[listsOfpolys//First,xvars,coeffs];
	STIIMapleScoreOrders[FirstScoredOrders,Drop[listsOfpolys,{1}],xvars,coeffs]
]


STIIMapleScoreOrders[ScoredOrders_,listsOfpolys_,xvars_,coeffs_]:=Module[{newScores},
	(*Print[ScoredOrders];*)
	(* New Scores *)
	newScores=STIIMapleScoreOrdersForPolys[ScoredOrders[[;;,1]],listsOfpolys//First,xvars,coeffs];
	(* Add Old Scores *)
	newScores=Table[nsc[[1]]-> nsc[[2]]+(nsc[[1]]/.ScoredOrders),{nsc,newScores}];
	(* Iterate *)
	STIIMapleScoreOrders[newScores,Drop[listsOfpolys,1],xvars,coeffs]
]


STIIMapleScoreOrders[ScoredOrders_,{},xvars_,coeffs_]:=ScoredOrders


(*STIIMapleScoreOrders[ScoredOrders_,listsOfpolys_,xvars_,coeffs_]:=Module[{orders,commonKeys},
]*)


(* Ask Maple to keep only LR orders and score them  *)
Clear[STIIMapleScoreOrdersForPolys];
STIIMapleScoreOrdersForPolys[orders_,polys_,xvars_,coeffs_]:=Module[{
polysAndPairs,
prefix,saveString,fileMapleCommand,
fileMapleOut,
rawResult
},

	
	polysAndPairs="["<>({polys,Subsets[polys,{2}]}//ToString[#,InputForm]&//stStripMmaContexts//StringTake[#,2;;-2]&)<>"]:";
	

	(* Main Maple Script *)
	prefix="with(MmaTranslator);\nread \""<>$SThyperIntPath<>"\";\n";
	prefix=prefix<>"vars:=["<>(ToString[xvars]//StringTake[#,2;;-2]&)<>"]:\n";

	prefix=prefix<>("unassign(`L`):L[{}]:="<>polysAndPairs<>"\n");
	prefix=prefix<>"cgReduction(L,"<>stStripMmaContexts[ToString[coeffs]]<>"):\n";
	
	prefix=prefix<>STIIMapleModCheckIntegrationOrder<>"\n"<>STIIMapleModFindAllOrders<>"\n"<>STIIMapleModFindOrdersWithScore<>"\n";
	prefix=prefix<>STIIMapleModScoreOrders<>"\n";
	prefix=prefix<>"orders:="<>(orders//ToString//StringReplace[#,{"{"->"[","}"->"]"}]&)<>":\n";
	
	
	
	prefix=prefix<>"finalRes:=ModScoreOrders(orders,L):\n";
	(*prefix//Print;*)
	
	(* Result String*)
	fileMapleOut=CreateFile[];
	saveString="fileRes:=\""<>fileMapleOut<>"\":\n";
	saveString=saveString<>"if FileTools[Exists](fileRes) then FileTools[Remove](fileRes): end if:\n FileTools[Text][Open](fileRes):\n";
	saveString=saveString<>"FileTools[Text][WriteString](fileRes,convert(finalRes,string)):\n";
	saveString=saveString<>"FileTools[Text][Close](fileRes):";
	
	(* Save Maple Command *)
	fileMapleCommand=CreateFile[];
	(*prefix<>"\n"<>saveString//Print;*)
	WriteString[fileMapleCommand,prefix<>"\n"<>saveString];
	Close[fileMapleCommand];
	
	Run[$MapleCommand<>" -q "<>fileMapleCommand];
	(*ReadString[fileMapleOut]//Print;*)
	ReadString[fileMapleOut]//StringReplace[#,{"\""->"",  "["->"{",  "]"-> "}","="-> "->" }]&//ToExpression

]


(* Ask Maple to compute all LR orders for a list of polynomials *)
Clear[STIIMapleFindAllLROrdersEndingWith];
STIIMapleFindAllLROrdersEndingWith[polys_,xvars_,coeffs_,xvarsend_]:=Module[{
polysAndPairs,
prefix,saveString,fileMapleCommand,
fileMapleOut,
rawResult,
ordersScores
},
	
	polysAndPairs="["<>({polys,Subsets[polys,{2}]}//ToString[#,InputForm]&//stStripMmaContexts//StringTake[#,2;;-2]&)<>"]:";

	(* Main Maple Script *)
	prefix="with(MmaTranslator);\nread \""<>$SThyperIntPath<>"\";\n";
	prefix=prefix<>"vars:=["<>(ToString[xvars]//StringTake[#,2;;-2]&)<>"]:\n";
	prefix=prefix<>"varsend:=["<>(ToString[xvarsend]//StringTake[#,2;;-2]&)<>"]:\n";
	

	prefix=prefix<>("unassign(`L`):L[{}]:="<>polysAndPairs<>"\n");
	prefix=prefix<>"cgReduction(L,"<>stStripMmaContexts[ToString[coeffs]]<>"):\n";
	(*prefix=prefix<>STIIMapleModCheckIntegrationOrder<>"\n"<>STIIMapleModFindAllOrders<>"\n"<>STIIMapleModFindOrdersWithScore<>"\n";*)
	prefix=prefix<>STIIMapleNestedPermutationsEndingWith<>"\n";
	
	
	prefix=prefix<>"costFun:=proc(Polys)\n numelems(Polys): end proc:"<>"\n";
	prefix=prefix<>"finalRes:=NestedPermutationsEndingWith(L,vars,varsend,costFun):\n";
	
	(*prefix//Print;*)
	
	(* Result String*)
	fileMapleOut=CreateFile[];
	saveString="fileRes:=\""<>fileMapleOut<>"\":\n";
	saveString=saveString<>"if FileTools[Exists](fileRes) then FileTools[Remove](fileRes): end if:\n FileTools[Text][Open](fileRes):\n";
	saveString=saveString<>"FileTools[Text][WriteString](fileRes,convert(finalRes,string)):\n";
	saveString=saveString<>"FileTools[Text][Close](fileRes):";
	
	savestring=prefix<>"\n"<>saveString;
	(*savestring//Print;*)
	
	(* Save Maple Command *)
	fileMapleCommand=CreateFile[];
	WriteString[fileMapleCommand,prefix<>"\n"<>saveString];
	Close[fileMapleCommand];
	
	(*fileMapleCommand//Print;*)
	(*{fileMapleCommand,fileMapleOut}//Print;*)
	
	Run[$MapleCommand<>" -q "<>fileMapleCommand];
	
	
	ordersScores=ReadString[fileMapleOut];
	If[ordersScores===EndOfFile,Return["error"],ordersScores=ordersScores//StringReplace[#,{"\""->"",  "["->"{",  "]"-> "}","="-> "->" }]&//ToExpression];
	Table[oS[[1]]->Total[oS[[2]]],{oS,ordersScores}]
]


(* Ask Maple to compute all LR orders for a list of polynomials *)
Clear[STIIMapleFindAllLROrders];
STIIMapleFindAllLROrders[polys_,xvars_,coeffs_]:=Module[{
polysAndPairs,
prefix,saveString,fileMapleCommand,
fileMapleOut,
rawResult,
ordersScores
},

	polysAndPairs="["<>({polys,Subsets[polys,{2}]}//ToString[#,InputForm]&//stStripMmaContexts//StringTake[#,2;;-2]&)<>"]:";

	(* Main Maple Script *)
	prefix="with(MmaTranslator);\nread \""<>$SThyperIntPath<>"\";\n";
	prefix=prefix<>"vars:=["<>(ToString[xvars]//StringTake[#,2;;-2]&)<>"]:\n";

	prefix=prefix<>("unassign(`L`):L[{}]:="<>polysAndPairs<>"\n");
	prefix=prefix<>"cgReduction(L,"<>stStripMmaContexts[ToString[coeffs]]<>"):\n";
	prefix=prefix<>STIIMapleModCheckIntegrationOrder<>"\n"<>STIIMapleModFindAllOrders<>"\n"<>STIIMapleModFindOrdersWithScore<>"\n";
	prefix=prefix<>STIIMapleNestedPermutations<>"\n";
	
	prefix=prefix<>"costFun:=proc(Polys)\n numelems(Polys): end proc:"<>"\n";
	prefix=prefix<>"finalRes:=NestedPermutations(L,vars,costFun):\n";
	prefix=prefix<>"finalRes1:=finalRes[..,1]:\n";
	prefix=prefix<>"finalRes1:=finalRes[..,2]:\n";
	(*prefix//Print;*)
	
	(* Result String*)
	fileMapleOut=CreateFile[];
	saveString="fileRes:=\""<>fileMapleOut<>"\":\n";
	saveString=saveString<>"if FileTools[Exists](fileRes) then FileTools[Remove](fileRes): end if:\n FileTools[Text][Open](fileRes):\n";
	saveString=saveString<>"FileTools[Text][WriteString](fileRes,convert(finalRes,string)):\n";
	saveString=saveString<>"FileTools[Text][Close](fileRes):";
	
	(* Save Maple Command *)
	fileMapleCommand=CreateFile[];
	WriteString[fileMapleCommand,prefix<>"\n"<>saveString];
	Close[fileMapleCommand];
	
	Run[$MapleCommand<>" -q "<>fileMapleCommand];
	(*ReadString[fileMapleOut]//Print;*)
	(*{fileMapleCommand,fileMapleOut}//Print;*)
	ordersScores=ReadString[fileMapleOut]//StringReplace[#,{"\""->"",  "["->"{",  "]"-> "}","="-> "->" }]&//ToExpression;
	Table[oS[[1]]->Total[oS[[2]]],{oS,ordersScores}]
]


(* ::Subsubsection::Closed:: *)
(*Additional LR functionalities*)


(* Ask Maple to suggest LR order *)
STIIMapleSuggestLROrder[polys_,coeffs_]:=Module[
{prefix,fileMapleOut,saveString,fileMapleCommand},

	(* Main Maple Script *)
	prefix="with(MmaTranslator);\nread \""<>$SThyperIntPath<>"\";\n";
	prefix=prefix<>("unassign(`L`):L[{}]:=["<>({polys,
	Subsets[polys,{2}]}//ToString[#,InputForm]&//stStripMmaContexts//StringTake[#,2;;-2]&)<>"]:\n");
	prefix=prefix<>"cgReduction(L,"<>stStripMmaContexts[ToString[coeffs]]<>"):\n";
	prefix=prefix<>"ord:=suggestIntegrationOrder(L):\n";
	prefix=prefix<>"finalRes:=ord:\n";
	(*prefix//Print;*)
	
	(* Result String*)
	fileMapleOut=CreateFile[];
	saveString="fileRes:=\""<>fileMapleOut<>"\":\n";
	saveString=saveString<>"if FileTools[Exists](fileRes) then FileTools[Remove](fileRes): end if:\n FileTools[Text][Open](fileRes):\n";
	saveString=saveString<>"FileTools[Text][WriteString](fileRes,convert(finalRes,string)):\n";
	saveString=saveString<>"FileTools[Text][Close](fileRes):";
	
	(* Save Maple Command *)
	fileMapleCommand=CreateFile[];
	WriteString[fileMapleCommand,prefix<>"\n"<>saveString];
	Close[fileMapleCommand];
	
	Run[$MapleCommand<>" -q "<>fileMapleCommand];
	(*ReadString[fileMapleOut]//Print;*)
	ReadString[fileMapleOut]//StringReplace[#,{"\""->"",  "["->"{",  "]"-> "}" }]&//ToExpression

]


(* Ask Maple to compute compatibility-graph reduction table *)
STIIMaplecgReduction[polys_,coeffs_]:=Module[{
prefix,saveString,fileMapleCommand,
fileMapleOut,
rawResult
},
	(* Maple Command *)
	prefix="with(MmaTranslator);\nread \""<>$SThyperIntPath<>"\";\n";
	prefix=prefix<>("unassign(`L`):L[{}]:=["<>({polys,
	Subsets[polys,{2}]}//ToString[#,InputForm]&//stStripMmaContexts//StringTake[#,2;;-2]&)<>"]:\n");
	prefix=prefix<>"cgReduction(L,"<>stStripMmaContexts[ToString[coeffs]]<>"):\n";
	prefix=prefix<>"finalRes:=L:\n";
	
	(* Result String*)
	fileMapleOut=CreateFile[];
	(* Print table on a fileMapleOut *)
	saveString=STIImaplePrintTable<>"\nPrintAsAssoc(\""<>fileMapleOut<>"\"):";
	
	
	fileMapleCommand=CreateFile[];
	WriteString[fileMapleCommand,prefix<>"\n"<>saveString];
	Close[fileMapleCommand];
	
	(* Launch Maple*)
	(*fileMapleCommand//Print;*)
	
	
	(*prefix<>"\n"<>saveString//Return;*)
	Run[$MapleCommand<>" -q "<>fileMapleCommand];
	ReadString[fileMapleOut]//StringReplace[#,{"\""->"",  "["->"{",  "]"-> "}" }]&//ToExpression

]


(* ::Subsubsection::Closed:: *)
(*Useful Maple Scripts*)


Clear[STIIMapleNestedPermutationsEndingWith];
STIIMapleNestedPermutationsEndingWith="NestedPermutationsEndingWith := proc(L::table, X, Tends::list, f::procedure)
    local Xlist, Xset, n, ans, count, extend, T, TendSets,
          S, x, Snew, candidates, k, perm, order, i;

    Xlist := convert(X, list);
    Xset := convert(Xlist, set);
    n := nops(Xlist);

    # Normalize endings as sets
    TendSets := [seq(convert(T, set), T in Tends)];

    # Check that all endings are contained in Xset
    for T in TendSets do
        if not T subset Xset then
            error \"An ending set is not contained in X\";
        end if;
    od;

    if not assigned(L[Xset]) then
        return [];
    end if;

    ans := table();
    count := 0;

    extend := proc(S::set, perm::list, order::list)
        local x, Snew, candidates, k, T, m, allowedEndings, C;

        k := nops(S);

        if k = 0 then
            count := count + 1;
            ans[count] := [perm, order];
            return;
        end if;

        candidates := {};

        # Union of the candidates allowed by each possible ending
        for T in TendSets do
            m := nops(T);

            if k > n - m then
                C := S intersect T;
            else
                C := S minus T;
            end if;

            candidates := candidates union C;
        od;

        for x in candidates do
            Snew := S minus {x};

            if Snew = {} then
                extend(
                    Snew,
                    [x, op(perm)],
                    order
                );
            elif assigned(L[Snew]) then
                extend(
                    Snew,
                    [x, op(perm)],
                    [f(L[Snew]), op(order)]
                );
            end if;
        od;
    end proc;

    extend(Xset, [], [f(L[Xset])]);

    return [seq(ans[i], i = 1 .. count)];
end proc:"


STIIMapleNestedPermutations="NestedPermutations := proc(L::table, X::list, f::procedure)
    local ans, n, extend, x, Snew;

    ans := [];
    n := nops(X);

    extend := proc(prefix::list, S::set, order::list)
        local x, Snew;

        if nops(prefix) = n then
            ans := [op(ans), [prefix, order]];
            return;
        end if;

        for x in X do
            if not member(x, S) then
                Snew := S union {x};

                if assigned(L[Snew]) then
                    extend(
                        [op(prefix), x],
                        Snew,
                        [op(order), f(L[Snew])]
                    );
                end if;
            end if;
        od;
    end proc;

    for x in X do
        if assigned(L[{x}]) then
            extend(
                [x],
                {x},
                [f(L[{x}])]
            );
        end if;
    od;

    return ans;
	(*return [indices(ans,'pairs')]:*)
end proc:"


(* Define a procedure to print a table in Mathematica Format on a file *)
STIImaplePrintTable="PrintAsAssoc := proc(fname)
    local f, k, keys, n, i, val;
    f := fopen(fname, WRITE);
    keys := [indices(L, 'nolist')];
    n := nops(keys);
    fprintf(f, \"<|\");
    for i to n do
        k := keys[i];
        val := L[k]; print(k):print(val):
        fprintf(f, \"\\\"%a\\\" -> %a\", k, val);
        if i < n then fprintf(f, \", \") end if;
    end do;
    fprintf(f, \"|>\");
    fclose(f);
end proc:"


(* Define a procedure to find LR orders with Scores *)


STIIMapleModFindOrdersWithScore="ModFindOrdersWithScore:= proc(i,vars)
local result,L:
      L:=table():L[{}]:=[convert(i[1],set),convert(i[2],set)]:
      cgReduction(L,{}):
      result:=ModCheckIntegrationOrder(L,vars):
      return result:
end proc:
"


(* Score a list of orders, keeping only those that are LR with respect to cgTable L  *)


STIIMapleModScoreOrders="ModScoreOrders := proc (orders,L::table, $)
local n, sub, goodOrders, ordering,res:
      goodOrders:=table():
      for ordering in orders do
      res:=ModCheckIntegrationOrder(L,convert(ordering,list)):
          if res[1] then goodOrders[convert(ordering,list)]:=res[2]: end if:
      end do:
      return [indices(goodOrders,'pairs')]:
end proc:
"


(* Given cgReduction table, find all orders and score them *)


STIIMapleModFindAllOrders="ModFindAllOrders := proc (L::table, vars::list, $)
local n, sub, goodOrders, ordering,res:
      goodOrders:=table():
      for ordering in Iterator:-Permute(vars) do
      res:=ModCheckIntegrationOrder(L,convert(ordering,list)):
          if res[1] then goodOrders[convert(ordering,list)]:=res[2]: end if:
      end do:
      return [indices(goodOrders,'pairs')]:
end proc:
"


(* Check if an order is LR and compute its score *)


STIIMapleModCheckIntegrationOrder="ModCheckIntegrationOrder := proc (L::table, vars::list, $)
local n, sub, score:
score:=0:
        for n from 1 to numelems(vars) do
		sub := {op(vars[1 .. n-1])}:
                if not assigned('L'[sub]) then
                        error \"No reduction data for %1 not available\", sub:
                elif max(map(degree, L[sub][1], vars[n]))>1 then
                        return [false,score]:
                end if:
                score:=score+numelems(L[sub][1]);
        end do:
        return [true,score]:
end proc:"


(* ::Section:: *)
(*SubTropicaII*)


(* ::Text:: *)
(*Pipeline:*)
(*1) We start with an Euler Integrand \mathcal{I} *)
(*2) We apply NP continuation -> combination of NP-continued integrands with the geometric property*)
(*3) We apply the subtraction formula -> combination of renormalized face integrands*)
(*4) Each renormalized integrand is a combination of counter-terms integrands, to be integrated in the same order*)


(* ::Subsection::Closed:: *)
(*Setup Kernels*)


STIISetupKernels[nkernels_Integer] := Module[{currentDir, kernelsAlreadyOk,
        effectiveNKernels = nkernels},
    currentDir = Directory[];

    If[TrueQ[$HyperIntroduceAlgebraicLetters] &&
       !TrueQ[$STFindRootsParallelSafe] &&
       effectiveNKernels > 1,
        effectiveNKernels = 1];

    (* Serial mode: close any subkernels, skip all distribution *)
    If[effectiveNKernels <= 1,
        If[Length[Kernels[]] > 0, CloseKernels[]];
        $STActiveKernelCount = 0;
        $STRequestedKernelCount = effectiveNKernels;

        (* Still need the job-tracking directory *)
        If[DirectoryQ[$STJobTrackingDir], DeleteDirectory[$STJobTrackingDir, DeleteContents -> True]];
        CreateDirectory[$STJobTrackingDir];
        Put[{}, $STCompletedJobsLog];

        $KernelSetupQ = True;
        Return[];
    ];

    $STRequestedKernelCount = effectiveNKernels;

    kernelsAlreadyOk = $KernelSetupQ &&
                       $STActiveKernelCount === effectiveNKernels &&
                       Length[Kernels[]] >= effectiveNKernels;

    If[!kernelsAlreadyOk,
        (* Fresh launch: close any stale kernels and start exactly effectiveNKernels *)
        CloseKernels[];
        LaunchKernels[effectiveNKernels];
        $STActiveKernelCount = effectiveNKernels;

        (* Sync working directory + search paths on new kernels.  The path sync
           must precede the parallel Get below so that SubTropica.wl's
           FF/SPQR auto-detect inside the subkernel can locate the packages. *)
        With[{cd = currentDir, mainPath = $Path, mainLibPath = $LibraryPath},
            ParallelEvaluate[
                SetDirectory[cd];
                $Path = mainPath;
                $LibraryPath = mainLibPath]];

        (* SubTropica` has \[Tilde]4.7 K symbols.  DistributeDefinitions["SubTropica`"]
           costs \[Tilde]20 s on a 13-subkernel pool because the symbol-table install on
           each subkernel scales with N_symbols * N_kernels regardless of how
           little data each definition holds.  Parsing SubTropica.wl locally on
           each subkernel takes \[Tilde]1 s in parallel when subkernel mode
           skips the dep-probe banner and persisted-config read (see the
           Global`$STSubkernelMode gate in SubTropica.wl).  HyperIntica` is
           loaded as a BeginPackage dependency by the same Get, so it comes
           up on subkernels in the same step.  Block[Print] suppresses the
           per-subkernel splash. *)
        With[{stPath = FileNameJoin[{$SubTropicaInstallDir, "SubTropica.wl"}]},
            ParallelEvaluate[
                Global`$STSubkernelMode = True;
                Block[{Print = (Null)&},
                    Off[General::shdw]; Off[Integrate::shdw];
                    Quiet[Get[stPath], {Integrate::shdw, General::shdw}]]]];

        (* Overlay any runtime HyperIntica state that the main kernel mutated
           between package load and STSetupKernel.  SubTropica's parallel Get
           above brought up a fresh HyperIntica context on each subkernel; this
           DistributeDefinitions imprints the main kernel's current values
           (algebraic letter counters, caches, ...) on top of that.  Cost
           is \[Tilde]0.4 s. *)
        DistributeDefinitions["HyperIntica`"];
        DistributeDefinitions["SubTropicaII`"];

        (* Mirror master-resolved $*Path / $*Command / $*PolymakeFraction
           values to every subkernel.  Required because subkernels skip
           SubTropica.wl's persisted-config read and dependency-probe banner
           in subkernel mode (Global`$STSubkernelMode = True), which means
           any user-side ConfigureSubTropica or runtime mutation that hadn't
           been persisted yet would otherwise be lost.  This is also
           important for paths that auto-discover differently between master
           and subkernel (e.g., a custom HyperFLINT install that the master
           knows about via persisted config but a fresh subkernel Get does
           not). *)
        With[{
              polymakeCmd  = $PolymakeCommand,
              ginshCmd     = $GinshCommand,
              mapleCmd     = $MapleCommand,
              hyperIntPath = $SThyperIntPath,
              hfPath       = $STHyperFlintPath,
              hfDataPath   = $STHyperFlintDataPath,
              pythonCmd    = $PythonCommand,
              polymakeFrac = $PolymakeConcurrencyFraction,
              fiestaPath   = $FIESTAPath,
              amflowPath   = $AMFlowPath,
              literedPath  = $LiteRedPath,
              liteibpPath  = $LiteIBPPath,
              firePath     = $FIREPath,
              feyntropPath = $FeyntropPath,
              ffPath       = $FiniteFlowPath,
              algLet       = $HyperIntroduceAlgebraicLetters,
              algSafe      = $STFindRootsParallelSafe,
              ffOn         = $UseFFPolynomialQuotient,
              checkDv      = $HyperInticaCheckDivergences},
            ParallelEvaluate[
                $PolymakeCommand                = polymakeCmd;
                $GinshCommand                   = ginshCmd;
                $MapleCommand                   = mapleCmd;
                $SThyperIntPath                 = hyperIntPath;
                $STHyperFlintPath               = hfPath;
                $STHyperFlintDataPath           = hfDataPath;
                $PythonCommand                  = pythonCmd;
                $PolymakeConcurrencyFraction    = polymakeFrac;
                $FIESTAPath                     = fiestaPath;
                $AMFlowPath                     = amflowPath;
                $LiteRedPath                    = literedPath;
                $LiteIBPPath                    = liteibpPath;
                $FIREPath                       = firePath;
                $FeyntropPath                   = feyntropPath;
                $FiniteFlowPath                 = ffPath;
                $HyperIntroduceAlgebraicLetters = algLet;
                $STFindRootsParallelSafe        = algSafe;
                $UseFFPolynomialQuotient        = ffOn;
                $HyperInticaCheckDivergences    = checkDv]];

        (* Suppress FrontEndObject::notavail on sub-kernels running headlessly *)
        ParallelEvaluate[Off[FrontEndObject::notavail]];
        (*  mirror the main-kernel Off[] calls on sub-kernels so that
           the $Failed[rays] cascade warnings are silenced there too *)
        ParallelEvaluate[Off[Table::iterb]; Off[Part::partd]];

        (* Static-definition broadcasts \[LongDash] cold-launch only.  These
           function/symbol values don't change after package load, so once
           the subkernel has them they stay valid for the lifetime of the
           pool.  Doing this on every warm STSetupKernel was \[Tilde]1 s of
           wasted DistributeDefinitions overhead per call (per
           profile_v1.1.6.3 round-2 + adversarial review). *)
        DistributeDefinitions[$STJobTrackingDir, $STCompletedJobsLog];
        DistributeDefinitions[STLaunchHyperInticaAllKernelIntegrator];
        DistributeDefinitions[STIILaunchHyperInticaAllKernelIntegrator];
        DistributeDefinitions[STLaunchHyperInticaAll];
        DistributeDefinitions[STIILaunchHyperInticaAll];

    ,
        (* Kernels already live with correct count; skip the expensive restart.
           Just re-sync the directory in case it changed.  No cache reset
           here \[LongDash] STSetupKernel is called multiple times inside a
           single STIntegrate run (gauge scoring, LR search, integration
           coordination), so clearing caches in the warm path would wipe
           state mid-integration.  Use the public STResetKernelCaches[]
           between STIntegrate calls when a clean slate is desired. *)
        ParallelEvaluate[SetDirectory[#]]&[currentDir];
    ];

    (* Always recreate the job-tracking directory so each integration run starts
       with a clean slate (required by STLaunchHyperIntica). *)
    If[DirectoryQ[$STJobTrackingDir], DeleteDirectory[$STJobTrackingDir, DeleteContents -> True]];
    CreateDirectory[$STJobTrackingDir];
    Put[{}, $STCompletedJobsLog];

    (* Propagate HyperIntica flags from main kernel to all sub-kernels *)
    With[{checkDiv = $HyperInticaCheckDivergences},
        ParallelEvaluate[$HyperInticaCheckDivergences = checkDiv]];

    $KernelSetupQ = True;
]


(* ::Subsection::Closed:: *)
(*Step-by-Step*)


(* ::Subsubsection::Closed:: *)
(*Basic Functions*)


(* Prepare data to run HyperIntica's polynomial reduction *)
(* Construct a list of the polynomial and pairs appearing in each counter-term of a locally finite integrand *)
(* The locally finite integrand can be passed directly as a Mathematica expression (not a list) *)

Clear[STIIpreparePolysAndPairs];
STIIpreparePolysAndPairs[totintegrand_,vars_]:=Module[{prefix,letters,den,toIntegrate,toInt,expanded,listsOflistsPolys,pairsOfPolys,factors},

	(* If the integrand is 0, nothing to do. This will be recognized later when finding linearly reducible orders *)
	If[totintegrand===0,Return[{}]];

	(* Expand into terms of the form R0 Log[R1]^n1 Log[R2]^n2 ... , with rational functions Ri *)
	expanded=totintegrand//Expand;
	(* Produce a list of the various terms *)
	If[Not[Head[expanded]===Plus],expanded={expanded},expanded=List@@expanded];
	(* Loop over individual terms, extract polynomial factors appearing in Ri (i>0) and in the denominator of R0 *)
	listsOflistsPolys=Table[
		(* Arguments of the logs *)
		letters=Cases[integrand,Log[a_]:> a,\[Infinity]]//DeleteDuplicates;
		(* Ignore constants *)
		letters//DeleteCases[ #,v_/;Count[v,Alternatives@@vars]===0]&;
		(* Denominator of prefactor *)
		den=Denominator[integrand/.Log[___]->1];
		(* Extract irreducible factors *)
		(**)
		If[letters==={},letters={{}}];
		factors=Join[letters,{den}]//STgetIrreducibleFacs//DeleteDuplicates//DeleteCases[ #,v_/;Count[{v},Alternatives@@vars,\[Infinity]]===0]&;
		factors
		,
		{integrand,expanded}
	]//DeleteDuplicates;
	(* Write pairs of polys for each integrand *)
	(*pairsOfPolys=Table[{polys,Subsets[polys,{2}]},{polys,listsOflistsPolys}];*)
	pairsOfPolys=Table[{polys,{}},{polys,listsOflistsPolys}];
	(* Find an order which works for all integrands, if none available find order that works for some of the first *)
	pairsOfPolys
]



(* ::Subsubsection:: *)
(*Setup directories (GS)*)


(* ::Text:: *)
(* Construct a directory to evaluate ONE subtraction formula at a given eps-order.*)
(**)
(* Optional arguments : *)
(*	 A string identifier for the directory*)
(*	construct polynomial and pairs for cg reduction?*)
(*	*)


STClearDirectories[id_:""]:=Module[{NPfilenames},
	NPfilenames=FileNames["integrands/"<>id<>"*"];
	DeleteDirectory[#,DeleteContents->True]&/@NPfilenames
]


STIIfastEpsSeriesCoefficient[integrand_,xvars_,order_]:=Module[{coefMonsPols,shortPolys,shortDefinitions,localPolyName,dummyIntegrand},
	coefMonsPols=STtoCoeffMonPols[integrand,xvars];
	shortPolys=Table[localPolyName[i],{i,1,coefMonsPols[[-1,1]]//Length}];
	
	shortDefinitions=Table[localPolyName[i]->coefMonsPols[[-1,1,i]],{i,1,coefMonsPols[[-1,1]]//Length}];
	
	dummyIntegrand=coefMonsPols[[1]](Times@@(xvars^coefMonsPols[[2]]))(Times@@(shortPolys^coefMonsPols[[-1,2]]));
	SeriesCoefficient[dummyIntegrand,{eps,0,order}]/.shortDefinitions
]


STIIfastEpsSeries[integrand_,xvars_,minOrder_,maxOrder_]:=Module[
{coefMonsPols,shortPolys,shortDefinitions,localPolyName,dummyIntegrand,
 series,seriesCoeffs,poleOrder,padFromSeries},
	(* Build a length-(maxOrder-minOrder+1) list whose k-th entry is the
	   \[CurlyEpsilon]^(minOrder + k - 1) coefficient of `s`.  SeriesData truncates trailing
	   zeros from high orders, and may start at poleOrder > minOrder when
	   the series has no pole that deep.  The original code used
	   PadLeft[seriesCoeffs, maxOrder-minOrder+1], which is correct only
	   when poleOrder == maxOrder; when poleOrder == minOrder (the
	   pure-pole counter-terms produced by tropical subtraction, e.g.
	   1/(2 eps^2) for a 2-loop triangle), PadLeft misplaces the pole
	   coefficient onto the \[CurlyEpsilon]^maxOrder slot and the poles disappear from
	   the final answer.  Fix: pad leading zeros for [minOrder, poleOrder-1]
	   and trailing zeros for [poleOrder+Length-1, maxOrder].  When
	   `s` is an \[CurlyEpsilon]-independent expression, Series[] returns a bare
	   value (not a SeriesData); we treat it as a length-1 series at
	   order 0. *)
	padFromSeries[s_]:=(
		If[Head[s]===SeriesData,
			{seriesCoeffs,poleOrder}={s[[3]],s[[4]]},
			{seriesCoeffs,poleOrder}={{s},0}];
		Join[
			ConstantArray[0,Max[0,poleOrder-minOrder]],
			PadRight[seriesCoeffs,Max[0,maxOrder-poleOrder+1]]]);

	If[Length[xvars]===0,
		(* Counter-term completely localised \[LongDash] no polynomial reconstruction
		   needed, Series operates directly on the integrand. *)
		padFromSeries[Series[integrand,{eps,0,maxOrder}]]
		,
		(* Non-localised integrand: use the short-name polynomial trick so
		   Series doesn't expand the (possibly large) base polynomials. *)
		coefMonsPols=STtoCoeffMonPols[integrand,xvars];
		shortPolys=Table[localPolyName[i],{i,1,Length[coefMonsPols[[-1,1]]]}];
		shortDefinitions=Table[localPolyName[i]->coefMonsPols[[-1,1,i]],
			{i,1,Length[coefMonsPols[[-1,1]]]}];
		dummyIntegrand=coefMonsPols[[1]]*
			(Times@@(xvars^coefMonsPols[[2]]))*
			(Times@@(shortPolys^coefMonsPols[[-1,2]]));
		padFromSeries[Series[dummyIntegrand,{eps,0,maxOrder}]]/.shortDefinitions
	]
]


Clear[STIIsetupDirectorySubtraction]
Options[STIIsetupDirectorySubtraction]={"MethodPolysAndPairs"->"Fast","LRordersQ"->False};
STIIsetupDirectorySubtraction[prefactor_,subtraction_,minOrder_,maxOrder_,coeffs_,id_:"subtraction",polysAndPairsQ_:True,OptionsPattern[]]:=Module[
{prefix,i,renintegrands,integrands,rnint,
variablesAndCoeffs,nfaces,polysAndPairs,directoryName,
dataFile,oldDirectory,messageFunction,
fromSetToList,exps,pols}
,

	messageFunction=If[$Notebooks,PrintTemporary,Print];

	(* How many faces *)
	nfaces=subtraction//Length;
	(* We treat the renormalized integrand of each face individually *)
	(*renintegrands= Table[SeriesCoefficient[face[[1]],{eps,0,order}],{face,subtraction}];*)
	renintegrands= Table[STIIfastEpsSeries[#,face[[2]],minOrder,maxOrder]&/@face[[1]],{face,subtraction}];

	(* integrands of each renormalized integrand *)
	integrands=Table[rnint,{rnint,renintegrands}];


	(* Prepare Face sub-Directories *)
	(* Loop over faces *)
	Do[
		Print["eps-order: "<>ToString[ord]];
		(* Polys and pairs of each ren-integrand *)
		If[OptionValue["MethodPolysAndPairs"]==="Fast",
			(*If[polysAndPairsQ,
					polysAndPairs=Table[
						exps=({(STtoCoeffMonPols[#,face[[2]]]//Last//Last),{}})&/@ face[[1]];
						exps=exps//Position[#,v_/;Count[{v},eps,\[Infinity]]>0 || v<0,1]&//Flatten;
						pols=({(STtoCoeffMonPols[#,face[[2]]]//Last//First),{}})&/@ face[[1]];
						pols[[exps]]
						,
						{face,subtraction}
					];
			],*)
			If[polysAndPairsQ,
					polysAndPairs=Table[
						Table[
							{pols,exps}=ct//(STtoCoeffMonPols[#,face[[2]]])&//Last;
							exps=exps//Position[#,v_/;Count[{v},eps,\[Infinity]]>0 || v<0 ,1]&//Flatten;
							{pols[[exps]],{}}
							,{ct,face[[1]]}
						]
						,{face,subtraction}
					];
			],
			If[polysAndPairsQ,
				polysAndPairs=Table[
					rnint=renintegrands[[i,;;,ord-minOrder+1]];
					rnint//Total//STIIpreparePolysAndPairs[#, subtraction[[i,2]] ]&
					,
					{i,1,subtraction//Length}
				];
			];
		];


		Do[

			(* Variables & coeffs of each renormalized-integrand *)
			variablesAndCoeffs = Table[{e[[2]],coeffs},{e,subtraction}];

			directoryName="integrands/"<>id<>"/ord_"<>ToString[ord]<>"_face_"<>ToString[i];
			If[Not[FileExistsQ[directoryName]],
				CreateDirectory[directoryName, CreateIntermediateDirectories -> True];
			];
			(* prefix pointing to the directory and face*)
			prefix=directoryName<>"/";
			(* Variables & Coefficients *)
			dataFile=prefix<>"vars.m";
			If[FileExistsQ[dataFile],DeleteFile[dataFile]];
			Put[variablesAndCoeffs[[i]],dataFile];
			(* Polys and pairs *)
			If[polysAndPairsQ,
				dataFile=prefix<>"polys.m";
				If[FileExistsQ[dataFile],DeleteFile[dataFile]];
				Put[polysAndPairs[[i]],dataFile];
			];
			(* Integrands *)
			dataFile=prefix<>"counter_terms_integrands.m";
			If[FileExistsQ[dataFile],DeleteFile[dataFile]];
			Put[integrands[[i,;;,ord-minOrder+1]],dataFile];
			(* eps-order (convenient later) *)
			dataFile=prefix<>"epsorder.m";
			If[FileExistsQ[dataFile],DeleteFile[dataFile]];
			Put[eps^ord,dataFile];
			(* prefactor (convenient later) *)
			dataFile=prefix<>"prefactor.m";
			If[FileExistsQ[dataFile],DeleteFile[dataFile]];
			Put[prefactor,dataFile];
			(* counter-terms results directory *)
			dataFile=prefix<>"partial_results";
			If[DirectoryQ[dataFile],DeleteDirectory[dataFile]];
			CreateDirectory[dataFile];
			(* LR created during expansion *)
			If[OptionValue["LRordersQ"],
				dataFile=prefix<>"bestOrder.m";
				If[FileExistsQ[dataFile],DeleteFile[dataFile]];
				Put[subtraction[[i,2]],dataFile];
			]
		,{i,1,nfaces}
		]
	,{ord,minOrder,maxOrder}
	];
]


Clear[STIIsetupDirectoryExpansion];
Options[STIIsetupDirectoryExpansion]={"MethodPolysAndPairs"->"Fast","LRordersQ"->False};
STIIsetupDirectoryExpansion[expansion_, order_, vars_, coeffs_ , id_:"NP", opts:OptionsPattern[] ]:=Module[
{
messageFunction,
prefactorEpsOrder,subtractionLeadingPoleOrder,minOrder, maxOrder,
results,
ord, cnt,
cntId,
cntIdlast,safeLoop,minOrderLast,maxOrderLast,fileLastResult,
dummy
},
	messageFunction=If[$Notebooks,Print,Print];

	(* Set up subtraction directory *)
	Print["---------------------------------"];
	Print["Setting up subtraction directories ! "];
	Do[
		cntId=id<>"_"<>ToString[cnt];
		messageFunction["cnt: "<>ToString[cnt]];
		prefactorEpsOrder=expansion[[cnt,1]]//STleadingLaurentOrder[#,eps]&;
		subtractionLeadingPoleOrder=-(Length[vars]-Min[Length/@expansion[[cnt,2,;;,2]]]);
		minOrder=subtractionLeadingPoleOrder;
		If[prefactorEpsOrder>0,
			maxOrder=order-prefactorEpsOrder,
			maxOrder=order-prefactorEpsOrder
		];
		STIIsetupDirectorySubtraction[expansion[[cnt,1]],expansion[[cnt,2]] , minOrder,maxOrder, coeffs , cntId, True, opts ];
		,
		{cnt,1,expansion//Length}
	];

]
(*STsetupDirectorySubtraction[subtraction_,order_,coeffs_,id_:"subtraction",polysAndPairsQ_:True]:=Module[*)


(* ::Subsubsection::Closed:: *)
(*Find orders *)


(* Enumerate all linearly reducible orders *)
Clear[STIIfindLinearlyReducibleOrders];
Options[STIIfindLinearlyReducibleOrders] = {
    "Strategy" -> "HighestEpsOrder",
    "Engine" -> "HyperIntica",
    Heuristic -> "LeafCountLinear",
    FindRoots -> False
};

STIIfindLinearlyReducibleOrders[id_:"NP", opts:OptionsPattern[]] := Module[{},
    Switch[
        OptionValue["Strategy"],
        "HighestEpsOrder", STIIfindLinearlyReducibleOrdersHighestEpsOrder[id,opts ]
    ]
];


Clear[STIIfindLinearlyReducibleOrdersHighestEpsOrder];
(* Find LR-orders in Highest eps-order, copy to lower ones *)
Options[STIIfindLinearlyReducibleOrdersHighestEpsOrder] = {
    Heuristic -> "LeafCountLinear",
    "Engine" -> "HyperIntica",
    FindRoots -> False
};

STIIfindLinearlyReducibleOrdersHighestEpsOrder[id_:"NP", OptionsPattern[]] := Module[
{
  NPfilenamesUpToEps, NPfilenamesEpsOrders,
  polysAndPairs, xvars, coeffs,linearlyReducibleOrders, bestOrder,
  failed = {},
  stringReplacerFunction,
  fileEpsMaxEpsOrder, fileEpsOrders, file,
  counter = 0, total,
  startTime
},

	NPfilenamesEpsOrders = STListDirectoriesNP[id];
	NPfilenamesUpToEps = NPfilenamesEpsOrders[[;; , 1]] // DeleteDuplicates;
	total = Length[NPfilenamesUpToEps];
    
	Print["---------------------------------"];
  
	startTime = AbsoluteTime[];
	Do[
		counter++;
		fileEpsOrders = Cases[NPfilenamesEpsOrders, {fileUpToEps, a_, b_} :> {a, b}] // SortBy[#, -#[[1]] &] &;
		file = fileEpsOrders // SortBy[#, (StringCases[#, "ord_" ~~ x___ ~~ "_face" :> ToExpression[x]]) &] & // Last;
		file = file[[2]];
		file // Print["Working on max eps-order : ", #] &;
		polysAndPairs = Get[file <> "/polys.m"];
		{xvars,coeffs} = Get[file <> "/vars.m"];      
		
		If[xvars === {} || polysAndPairs === {},
			Print["--- No integration required!"];
			linearlyReducibleOrders = "no_integration_required";
			bestOrder = "no_integration_required";
			,
			Module[{espResult, rootPolys = {}},
				espResult = Switch[OptionValue["Engine"],
					"HyperIntica",
					STEspressoFubini[Join[#, xvars] & /@ (polysAndPairs[[;; , 1]]), xvars, Heuristic -> OptionValue[Heuristic], FindRoots -> OptionValue[FindRoots]],
					"HyperInt",
					Print["Using Maple!"];
					{(polysAndPairs[[;; , 1]])//STIIMapleFindLROrders[#,xvars,coeffs]&//SortBy[#,Last]&//First//First,0 (* dummy score *)}
					
				];
				
				
				If[OptionValue[FindRoots],
					{bestOrder, score} = espResult[[1]];
					rootPolys = espResult[[2]];
					,
					{bestOrder, score} = espResult;
				];
				If[bestOrder === NOLR||bestOrder==={},
					failed = Join[failed, {file}];
					Message[STEspressoFubini::noorder, file, xvars];
					Abort[];
				];
			
				(* deg-3 NOLR detection only; HyperInt introduces Wm[i]/Wp[i] on demand via $HyperIntroduceAlgebraicLetters. *)
				If[OptionValue[FindRoots] && rootPolys =!= {},
					STApplyRootFactoring[polysAndPairs, rootPolys, xvars, bestOrder];
					Print["  [FindRoots] Identified ", Length[DeleteDuplicates[rootPolys]]
					,
					" root polynomial(s) (HyperInt will introduce Wm/Wp letters)"];
				];
			];
		];
		
		Do[
			Put[bestOrder, file2 <> "/bestOrder.m"];
			,
			{file2, fileEpsOrders[[;; , 2]]}
		];
		,
		{fileUpToEps, NPfilenamesUpToEps}
	]
]


STEspressoFubini::noorder = "No linearly reducible integration order found for `1`. Either no such order exists, or it does not exist for the current choice of variables: `2`.";


(* ::Subsubsection::Closed:: *)
(*Launch HyperIntica *)


(* Enumerate all linearly reducible orders *)
Clear[STIILaunchHyperIntica];

Options[STIILaunchHyperIntica] = {
    "Strategy" -> "All",
    "LevelParallelism" -> "Face",
    "NumericsQ" -> False,
    "ShowIntegrands" -> False,
    "ClearCachesPerIntegrand" -> False,
    "ClearCachesMemoryThreshold" -> Infinity,
    "SelectFaces" -> All,
    "ReuseExistingResults" -> True,
    "UIComms" -> None,
    "Integrator" -> "HyperIntica"
};

STIILaunchHyperIntica::nokernel = "Parallelization -> \"PartialIntegrands\" requires subkernels, but KernelsAvailable was set to 0 or 1 (serial mode). Use Parallelization -> All or increase KernelsAvailable.";

STIILaunchHyperIntica[id_:"NP", OptionsPattern[]] := Module[
{NPfilenames, polysAndPairs, xvars, linearlyReducibleOrders,
successQ, failed = {}},
    If[OptionValue["Strategy"] === "PartialIntegrands" && Length[Kernels[]] === 0,
        Message[STLaunchHyperIntica::nokernel]; Abort[];
    ];
    Switch[OptionValue["Strategy"],
    All,
        failed = STLaunchHyperInticaAll[id,
            "NumericsQ"                  -> OptionValue["NumericsQ"],
            "ShowIntegrands"             -> OptionValue["ShowIntegrands"],
            "ClearCachesPerIntegrand"    -> OptionValue["ClearCachesPerIntegrand"],
            "ClearCachesMemoryThreshold" -> OptionValue["ClearCachesMemoryThreshold"],
            "SelectFaces"                -> OptionValue["SelectFaces"],
            "ReuseExistingResults"       -> OptionValue["ReuseExistingResults"],
            "UIComms"                    -> OptionValue["UIComms"],
            "Integrator"                 -> OptionValue["Integrator"]
        ]
    ];
    failed
]


$KernelSetupQ = False;
$STActiveKernelCount = 0;          (*  track running kernel count for reuse check *)
$STTropicalDataCache = <||>;       (*  in-memory cache for polymake tropical fan results *)
$PolymakeConcurrencyFraction = 0.5; (*  default fraction of CPU cores for concurrent polymake jobs *)
$STPolymakeProcess = None;          (*  handle to the currently running polymake batch bash process, for abort cleanup *)

(* $STFindRootsParallelSafe: when True (default as of v1.0.402),
   subkernels run FindRoots->True integrations with Block-scoped
   $HyperAlgebraicLetterCounter / Table, and the main-kernel aggregator
   merges the per-subkernel letter tables before SimplifyWithVieta /
   CanonicalizeAlgebraicLetters run. The scoring stage
   (STfindLinearlyReducibleOrders2 via STFasterFubini2) is already
   safe by construction: it uses a local rootCounter and does NOT
   mutate $HyperAlgebraicLetterCounter/Table, and STApplyRootFactoring
   is skipped when "ScanGauges" -> True. Integration is made safe by
   the JobIndex*$STFindRootsJobStride Block-scope inside
   STLaunchHyperInticaAllKernelIntegrator.

   Set False as an opt-out to restore the legacy force-serial path in
   STSetupKernel (closes all subkernels whenever
   $HyperIntroduceAlgebraicLetters is True).

   Only the "All" strategy is parallel-safe; BruteForce and
   PartialIntegrands have their own inline force-serial guards. *)
$STFindRootsParallelSafe = True;

(* Per-subkernel index stride for FindRoots letter allocation. Each
   parallel job gets a disjoint integer range [jobIndex*stride,
   (jobIndex+1)*stride) for Wm[i]/Wp[i] indices, guaranteeing that
   different subkernels' letters never collide when their tables are
   unioned on the main kernel. 10^6 is overkill (worst library topology
   uses ~10 letters) but costs nothing and leaves plenty of headroom. *)
$STFindRootsJobStride = 10^6;

(*  suppress harmless ANNOYING iterator/part warnings that fire when polymake
   returns $Failed (e.g. Table[...,{ray,$Failed[rays]}] and ray[[1]] on $Failed).
   General::stop is a follow-on meta-message that disappears automatically. *)
Off[Table::iterb];
Off[Part::partd];
$STJobTrackingDir = FileNameJoin[{$TemporaryDirectory, "STJobTracking"}];
$STCompletedJobsLog = FileNameJoin[{$TemporaryDirectory, "STJobTracking", "completed_jobs.m"}];


(* STResetKernelCaches[] mirrors what a kernel restart does for the
   caches that affect result correctness or performance, on both the
   master kernel and every active subkernel.  Called automatically inside
   STSetupKernel's warm path so re-using an already-launched pool is
   indistinguishable from a fresh launch.  The Long benchmark suite (22
   diagrams, back-to-back, single kernel session) has stable result
   hashes across runs, which is what validates that this set of clears
   is sufficient. *)
STResetKernelCaches[] := (
    Quiet[
        If[ValueQ[HyperIntica`$LinearFactorsCache],          HyperIntica`$LinearFactorsCache = <||>];
        If[ValueQ[$NoAlgebraicRootsContributions],           $NoAlgebraicRootsContributions = <||>];
        If[ValueQ[HyperIntica`$HyperAlgebraicLetterTable],   HyperIntica`$HyperAlgebraicLetterTable = <||>];
        If[ValueQ[HyperIntica`$HyperAlgebraicLetterCounter], HyperIntica`$HyperAlgebraicLetterCounter = 0];
        $STDispatchHFCount = 0; $STDispatchHICount = 0;
        $STDispatchHFTime  = 0.; $STDispatchHITime  = 0.;
        $STDispatchProfileLog = {};
        $STHyperFlintCallCount = 0; $STHyperFlintTotalTime = 0.;
        $STSetupDirCallCount = 0; $STSetupDirPutTime = 0.;
        ClearSystemCache[];
        Quiet[ForgetAllMemo[]];
        (* Clear the proportional-polynomial memos.  These accumulate
           inside the gauge-scoring ParallelTable and are NOT cleared by
           ForgetAllMemo (which only resets HyperIntica memos).  When an
           inner TimeConstrained interrupt fires (e.g. on an entry that
           overruns the 180s gate), the in-flight LR scan exits before
           reaching the natural ForgetProportionalPolynomialsQ cleanup
           at line 16251, leaving up to ~48 MB / 47k entries on master
           and ~237 MB across 13 subkernels.  Subsequent calls then see
           a 5-7x slowdown in the gauge phase.  ForgetProportionalPolynomialsQ
           and its LR sibling already do the parallel cleanup themselves;
           we just need to invoke them from the central reset. *)
        Quiet[ForgetProportionalPolynomialsQ[]];
        Quiet[ForgetProportionalPolynomialsQLR[]]];
    If[Length[Kernels[]] > 0,
        Quiet @ ParallelEvaluate[
            If[ValueQ[HyperIntica`$LinearFactorsCache],          HyperIntica`$LinearFactorsCache = <||>];
            If[ValueQ[$NoAlgebraicRootsContributions],           $NoAlgebraicRootsContributions = <||>];
            If[ValueQ[HyperIntica`$HyperAlgebraicLetterTable],   HyperIntica`$HyperAlgebraicLetterTable = <||>];
            If[ValueQ[HyperIntica`$HyperAlgebraicLetterCounter], HyperIntica`$HyperAlgebraicLetterCounter = 0];
            (* Mirror the master-side dispatch counter resets so
               aggregated counts after STIntegrate reflect only the
               current run, not accumulated history. *)
            $STDispatchHFCount = 0; $STDispatchHICount = 0;
            $STDispatchHFTime  = 0.; $STDispatchHITime  = 0.;
            $STDispatchProfileLog = {};
            $STHyperFlintCallCount = 0; $STHyperFlintTotalTime = 0.;
            $STSetupDirCallCount = 0; $STSetupDirPutTime = 0.;
            ClearSystemCache[];
            Quiet[ForgetAllMemo[]];
            $HistoryLength = 0;
            Share[]]];
    Share[];);

(* Eager-pool launch task placeholder (set at the end of SubTropica.wl
   when $STEagerKernelPool is True).  STSetupKernel's entry wrapper waits
   on it before doing any work so the first user-triggered call benefits
   from the already-warm pool.  See the kickoff block at the bottom of
   this file. *)
$STEagerLaunchTask = None;
$STEagerKernelPool = True;

(*  STSetupKernel skips CloseKernels+LaunchKernels when the correct
   number of kernels are already alive, saving ~5.5 s per redundant call,
   and clears caches on every reuse so warm-path semantics match a fresh
   launch.  The job-tracking directory is always recreated.

   nkernels <= 1 means fully serial mode: all subkernels are closed and
   everything runs on the main kernel.  ParallelTable/ParallelMap/etc.
   automatically fall back to their serial equivalents when no subkernels
   are alive.

   The public `STSetupKernel` is a thin wrapper that synchronizes with
   any pending eager-pool task before delegating to `stSetupKernelImpl`.
   The eager task itself calls `stSetupKernelImpl` directly so it does
   not deadlock on its own pending-task slot. *)
STSetupKernel[nkernels_Integer:3] := Module[{task},
    If[Head[$STEagerLaunchTask] === TaskObject,
        task = $STEagerLaunchTask;
        $STEagerLaunchTask = None;  (* clear before waiting so the wrapper
                                        doesn't recurse on its own task *)
        Quiet @ TaskWait[task];
        (* If the eager task already launched the requested pool size, the
           kernels are clean and ready \[LongDash] no need to fall through to
           stSetupKernelImpl, which would re-traverse the warm path and run
           a redundant STResetKernelCaches.  Fall through only when the
           eager run does not match (different size requested, or the
           background task aborted). *)
        If[TrueQ[$KernelSetupQ] &&
            $STActiveKernelCount === nkernels &&
            Length[Kernels[]] >= nkernels,
            Return[Null]]];
    stSetupKernelImpl[nkernels]]

stSetupKernelImpl[nkernels_Integer] := Module[{currentDir, kernelsAlreadyOk,
        effectiveNKernels = nkernels},
    currentDir = Directory[];

    (* FindRoots mode: LinearFactors mutates $HyperAlgebraicLetterCounter /
       $HyperAlgebraicLetterTable. In the legacy path these are not shared
       across subkernels, so parallel integration would lose letters or
       double-count indices \[LongDash] force serial. When
       $STFindRootsParallelSafe is True, the "All" strategy's per-kernel
       integrator (STLaunchHyperInticaAllKernelIntegrator) Block-scopes
       the counter/table with a unique JobIndex*$STFindRootsJobStride
       offset, and the main-kernel aggregator merges the tables before
       post-processing; in that mode the force is off. BruteForce and
       PartialIntegrands retain inline force-serial guards (see
       STLaunchHyperInticaBruteForce / ...PartialIntegrands). *)
    (* Silent force: FindRoots + flag-off falls back to serial by design
       now that FindRoots defaults to True. Users who want parallel
       algebraic letters opt in via $STFindRootsParallelSafe = True. *)
    If[TrueQ[$HyperIntroduceAlgebraicLetters] &&
       !TrueQ[$STFindRootsParallelSafe] &&
       effectiveNKernels > 1,
        effectiveNKernels = 1];

    (* Serial mode: close any subkernels, skip all distribution *)
    If[effectiveNKernels <= 1,
        If[Length[Kernels[]] > 0, CloseKernels[]];
        $STActiveKernelCount = 0;
        $STRequestedKernelCount = effectiveNKernels;

        (* Still need the job-tracking directory *)
        If[DirectoryQ[$STJobTrackingDir], DeleteDirectory[$STJobTrackingDir, DeleteContents -> True]];
        CreateDirectory[$STJobTrackingDir];
        Put[{}, $STCompletedJobsLog];

        $KernelSetupQ = True;
        Return[];
    ];

    $STRequestedKernelCount = effectiveNKernels;

    kernelsAlreadyOk = $KernelSetupQ &&
                       $STActiveKernelCount === effectiveNKernels &&
                       Length[Kernels[]] >= effectiveNKernels;

    If[!kernelsAlreadyOk,
        (* Fresh launch: close any stale kernels and start exactly effectiveNKernels *)
        CloseKernels[];
        LaunchKernels[effectiveNKernels];
        $STActiveKernelCount = effectiveNKernels;

        (* Sync working directory + search paths on new kernels.  The path sync
           must precede the parallel Get below so that SubTropica.wl's
           FF/SPQR auto-detect inside the subkernel can locate the packages. *)
        With[{cd = currentDir, mainPath = $Path, mainLibPath = $LibraryPath},
            ParallelEvaluate[
                SetDirectory[cd];
                $Path = mainPath;
                $LibraryPath = mainLibPath]];

        (* SubTropica` has \[Tilde]4.7 K symbols.  DistributeDefinitions["SubTropica`"]
           costs \[Tilde]20 s on a 13-subkernel pool because the symbol-table install on
           each subkernel scales with N_symbols * N_kernels regardless of how
           little data each definition holds.  Parsing SubTropica.wl locally on
           each subkernel takes \[Tilde]1 s in parallel when subkernel mode
           skips the dep-probe banner and persisted-config read (see the
           Global`$STSubkernelMode gate in SubTropica.wl).  HyperIntica` is
           loaded as a BeginPackage dependency by the same Get, so it comes
           up on subkernels in the same step.  Block[Print] suppresses the
           per-subkernel splash. *)
        With[{stPath = FileNameJoin[{$SubTropicaInstallDir, "SubTropica.wl"}]},
            ParallelEvaluate[
                Global`$STSubkernelMode = True;
                Block[{Print = (Null)&},
                    Off[General::shdw]; Off[Integrate::shdw];
                    Quiet[Get[stPath], {Integrate::shdw, General::shdw}]]]];

        (* Overlay any runtime HyperIntica state that the main kernel mutated
           between package load and STSetupKernel.  SubTropica's parallel Get
           above brought up a fresh HyperIntica context on each subkernel; this
           DistributeDefinitions imprints the main kernel's current values
           (algebraic letter counters, caches, ...) on top of that.  Cost
           is \[Tilde]0.4 s. *)
        DistributeDefinitions["HyperIntica`"];

        (* Mirror master-resolved $*Path / $*Command / $*PolymakeFraction
           values to every subkernel.  Required because subkernels skip
           SubTropica.wl's persisted-config read and dependency-probe banner
           in subkernel mode (Global`$STSubkernelMode = True), which means
           any user-side ConfigureSubTropica or runtime mutation that hadn't
           been persisted yet would otherwise be lost.  This is also
           important for paths that auto-discover differently between master
           and subkernel (e.g., a custom HyperFLINT install that the master
           knows about via persisted config but a fresh subkernel Get does
           not). *)
        With[{
              polymakeCmd  = $PolymakeCommand,
              ginshCmd     = $GinshCommand,
              mapleCmd     = $MapleCommand,
              hyperIntPath = $SThyperIntPath,
              hfPath       = $STHyperFlintPath,
              hfDataPath   = $STHyperFlintDataPath,
              pythonCmd    = $PythonCommand,
              polymakeFrac = $PolymakeConcurrencyFraction,
              fiestaPath   = $FIESTAPath,
              amflowPath   = $AMFlowPath,
              literedPath  = $LiteRedPath,
              liteibpPath  = $LiteIBPPath,
              firePath     = $FIREPath,
              feyntropPath = $FeyntropPath,
              ffPath       = $FiniteFlowPath,
              algLet       = $HyperIntroduceAlgebraicLetters,
              algSafe      = $STFindRootsParallelSafe,
              ffOn         = $UseFFPolynomialQuotient,
              checkDv      = $HyperInticaCheckDivergences},
            ParallelEvaluate[
                $PolymakeCommand                = polymakeCmd;
                $GinshCommand                   = ginshCmd;
                $MapleCommand                   = mapleCmd;
                $SThyperIntPath                 = hyperIntPath;
                $STHyperFlintPath               = hfPath;
                $STHyperFlintDataPath           = hfDataPath;
                $PythonCommand                  = pythonCmd;
                $PolymakeConcurrencyFraction    = polymakeFrac;
                $FIESTAPath                     = fiestaPath;
                $AMFlowPath                     = amflowPath;
                $LiteRedPath                    = literedPath;
                $LiteIBPPath                    = liteibpPath;
                $FIREPath                       = firePath;
                $FeyntropPath                   = feyntropPath;
                $FiniteFlowPath                 = ffPath;
                $HyperIntroduceAlgebraicLetters = algLet;
                $STFindRootsParallelSafe        = algSafe;
                $UseFFPolynomialQuotient        = ffOn;
                $HyperInticaCheckDivergences    = checkDv]];

        (* Suppress FrontEndObject::notavail on sub-kernels running headlessly *)
        ParallelEvaluate[Off[FrontEndObject::notavail]];
        (*  mirror the main-kernel Off[] calls on sub-kernels so that
           the $Failed[rays] cascade warnings are silenced there too *)
        ParallelEvaluate[Off[Table::iterb]; Off[Part::partd]];

        (* Static-definition broadcasts \[LongDash] cold-launch only.  These
           function/symbol values don't change after package load, so once
           the subkernel has them they stay valid for the lifetime of the
           pool.  Doing this on every warm STSetupKernel was \[Tilde]1 s of
           wasted DistributeDefinitions overhead per call (per
           profile_v1.1.6.3 round-2 + adversarial review). *)
        DistributeDefinitions[$STJobTrackingDir, $STCompletedJobsLog];
        DistributeDefinitions[STLaunchHyperInticaAllKernelIntegrator];
        DistributeDefinitions[STLaunchHyperInticaAll];

    ,
        (* Kernels already live with correct count; skip the expensive restart.
           Just re-sync the directory in case it changed.  No cache reset
           here \[LongDash] STSetupKernel is called multiple times inside a
           single STIntegrate run (gauge scoring, LR search, integration
           coordination), so clearing caches in the warm path would wipe
           state mid-integration.  Use the public STResetKernelCaches[]
           between STIntegrate calls when a clean slate is desired. *)
        ParallelEvaluate[SetDirectory[#]]&[currentDir];
    ];

    (* Always recreate the job-tracking directory so each integration run starts
       with a clean slate (required by STLaunchHyperIntica). *)
    If[DirectoryQ[$STJobTrackingDir], DeleteDirectory[$STJobTrackingDir, DeleteContents -> True]];
    CreateDirectory[$STJobTrackingDir];
    Put[{}, $STCompletedJobsLog];

    (* Propagate HyperIntica flags from main kernel to all sub-kernels *)
    With[{checkDiv = $HyperInticaCheckDivergences},
        ParallelEvaluate[$HyperInticaCheckDivergences = checkDiv]];

    $KernelSetupQ = True;
]


(* ::Subsubsection::Closed:: *)
(*Strategy: All*)


(* Snapshot of all HyperIntica session-state globals that affect LinearFactors /
   HyperInt / downstream integrator behavior. Built on the main kernel at
   dispatch time and passed through to each subkernel via the
   STLaunchHyperInticaAllKernelIntegrator "HyperSnapshot" option; the
   subkernel's Block rebinds every key explicitly.  This closes the
   class of parallel-safety bug where a subkernel's default value of a
   HyperIntica flag differs from the main kernel's current setting
   (v1.0.402 $STFindRootsParallelSafe regression: subkernels had
   $HyperIntroduceAlgebraicLetters = False, main had True, so
   LinearFactors silently zeroed algebraic contributions). Adding a new
   HyperIntica global that affects integration means adding a key here
   and to stApplyHyperSnapshotToBlock below.  Keys drop the leading "$". *)
stBuildHyperSnapshot[] := <|
    "HyperIntroduceAlgebraicLetters"  -> TrueQ[HyperIntica`$HyperIntroduceAlgebraicLetters],
    "NoAlgebraicRootsContributions"   -> TrueQ[HyperIntica`$NoAlgebraicRootsContributions],
    "HyperAlgebraicRoots"             -> TrueQ[HyperIntica`$HyperAlgebraicRoots],
    "HyperIgnoreNonlinearPolynomials" -> TrueQ[HyperIntica`$HyperIgnoreNonlinearPolynomials],
    "HyperWarnZeroed"                 -> TrueQ[HyperIntica`$HyperWarnZeroed],
    "HyperSplittingField"             -> HyperIntica`$HyperSplittingField,
    "HyperVerbosity"                  -> HyperIntica`$HyperVerbosity,
    "HyperInticaCheckDivergences"     -> TrueQ[HyperIntica`$HyperInticaCheckDivergences],
    "HyperInticaAbortOnDivergence"    -> TrueQ[HyperIntica`$HyperInticaAbortOnDivergence],
    "HyperEvaluatePeriods"            -> TrueQ[HyperIntica`$HyperEvaluatePeriods]
|>;


Options[STLaunchHyperInticaAllKernelIntegrator] = {
    "NumericsQ" -> False,
    "ShowIntegrands" -> False,
    "ClearCachesPerIntegrand" -> False,
    "ClearCachesMemoryThreshold" -> Infinity,
    "Integrator" -> "HyperIntica",
    "JobIndex" -> 0,      (*  0-based unique index across the ParallelMap dispatch.
                              Used to shard FindRoots letter-index space per subkernel. *)
    "JobStride" -> 0,     (*  Size of each subkernel's letter-index namespace. 0 means
                              "don't shard" (legacy behaviour, caller must ensure serial). *)
    "HyperSnapshot" -> <||>  (* Main-kernel HyperIntica session-state snapshot (see
                                stBuildHyperSnapshot); applied via Block inside the
                                subkernel so LinearFactors/HyperInt see the same flag
                                values the main kernel has.  Default <||> means "leave
                                subkernel's current values unchanged" (used for legacy
                                direct calls; parallel dispatch always passes a real
                                snapshot). *)
};

STLaunchHyperInticaAllKernelIntegrator[{faceDirectory_, ctId_, ctIntegrand_, LRorder_}, OptionsPattern[]] := Module[
{result, localLetterTable = <||>, successQ, rules, showIntegrands, startTime, endTime,
 duration, jobFile, jobInfo, logFile, memThreshold, integrator,
 jobIndex, jobStride, letterBase, runIntegrator, hyperSnap, snapVal},
    showIntegrands = OptionValue["ShowIntegrands"];
    memThreshold   = OptionValue["ClearCachesMemoryThreshold"];
    jobIndex       = OptionValue["JobIndex"];
    jobStride      = OptionValue["JobStride"];
    hyperSnap      = OptionValue["HyperSnapshot"];
    letterBase     = If[IntegerQ[jobIndex] && IntegerQ[jobStride] && jobStride > 0,
                        jobIndex * jobStride, 0];
    (* Resolve a snapshot key: use the snapshot's value if present, else the
       subkernel's current global. Keeps the Block below valid even when
       callers pass an empty snapshot (legacy / serial path). *)
    If[!AssociationQ[hyperSnap], hyperSnap = <||>];
    snapVal[key_String, default_] := Lookup[hyperSnap, key, default];
    startTime = AbsoluteTime[];

    Print["Working on : ", faceDirectory, " counter term id: ", ctId];

    (* Write job info to a file *)
    If[showIntegrands,
        jobFile = FileNameJoin[{$STJobTrackingDir, "kernel_" <> ToString[$KernelID] <> ".m"}];
        jobInfo = <|
            "Integrand" -> ctIntegrand,
            "Order" -> LRorder,
            "Face" -> faceDirectory,
            "CtId" -> ctId,
            "StartTime" -> startTime
        |>;
        Put[jobInfo, jobFile];
    ];

    (rules=OptionValue["NumericsQ"]);

    integrator = Switch[OptionValue["Integrator"],
        "HyperInt",   SThyperIntMaple,
        "HyperFLINT", stTimedHyperFlint,
        _,            HyperInt];

    (* runIntegrator computes the result in the current dynamic scope.  We
       Block-scope the HyperIntica letter-allocation state around this call
       so that when $HyperIntroduceAlgebraicLetters is True, this
       subkernel's Wm[i]/Wp[i] indices live in a disjoint range starting at
       letterBase, and the populated table is captured into
       localLetterTable for the main-kernel aggregator to union in.  When
       FindRoots is off this is a cheap no-op: Block rebinds the variables,
       but LinearFactors never touches them. *)
    runIntegrator[] := (
        If[LRorder==="no_integration_required",
    		If[rules===False,result =ctIntegrand,result=ctIntegrand/.rules],
    		If[rules===False,
    			result=integrator[ctIntegrand,LRorder]
    			,
    			(*result=HyperInt[ctIntegrand/. Dispatch[rules],LRorder]*)
    			result=integrator[Evaluate[ctIntegrand/. rules],LRorder]
    		];
    	];
    );

    (* Block-scope every HyperIntica flag the caller captured plus the
       per-subkernel letter counter/table.  Without this, a subkernel's
       default $HyperIntroduceAlgebraicLetters (False) would override the
       main kernel's current setting (True when FindRoots is active), and
       LinearFactors would silently zero deg-2 factors on each subkernel
       (v1.0.402 regression, bisected on the 2-mass alternating box). *)
    Block[{
        HyperIntica`$HyperAlgebraicLetterCounter = letterBase,
        HyperIntica`$HyperAlgebraicLetterTable   = <||>,
        HyperIntica`$HyperIntroduceAlgebraicLetters =
            snapVal["HyperIntroduceAlgebraicLetters",
                    HyperIntica`$HyperIntroduceAlgebraicLetters],
        HyperIntica`$NoAlgebraicRootsContributions =
            snapVal["NoAlgebraicRootsContributions",
                    HyperIntica`$NoAlgebraicRootsContributions],
        HyperIntica`$HyperAlgebraicRoots =
            snapVal["HyperAlgebraicRoots",
                    HyperIntica`$HyperAlgebraicRoots],
        HyperIntica`$HyperIgnoreNonlinearPolynomials =
            snapVal["HyperIgnoreNonlinearPolynomials",
                    HyperIntica`$HyperIgnoreNonlinearPolynomials],
        HyperIntica`$HyperWarnZeroed =
            snapVal["HyperWarnZeroed",
                    HyperIntica`$HyperWarnZeroed],
        HyperIntica`$HyperSplittingField =
            snapVal["HyperSplittingField",
                    HyperIntica`$HyperSplittingField],
        HyperIntica`$HyperVerbosity =
            snapVal["HyperVerbosity",
                    HyperIntica`$HyperVerbosity],
        HyperIntica`$HyperInticaCheckDivergences =
            snapVal["HyperInticaCheckDivergences",
                    HyperIntica`$HyperInticaCheckDivergences],
        HyperIntica`$HyperInticaAbortOnDivergence =
            snapVal["HyperInticaAbortOnDivergence",
                    HyperIntica`$HyperInticaAbortOnDivergence],
        HyperIntica`$HyperEvaluatePeriods =
            snapVal["HyperEvaluatePeriods",
                    HyperIntica`$HyperEvaluatePeriods]
    },
        runIntegrator[];
        localLetterTable = HyperIntica`$HyperAlgebraicLetterTable;
    ];

    (* successQ is `True` unless the integrator explicitly returned a
       `$Failed` sentinel (e.g. SThyperIntMaple returning $Failed when
       Maple's `hyperInt` emitted an un-roundtrippable `Root(...)`
       algebraic-root expression).  The aggregator at the face level
       checks for this and propagates failure cleanly instead of
       silently including `$Failed` in the symbolic Laurent series. *)
    successQ = !(result === $Failed || !FreeQ[result, $Failed]);
    endTime = AbsoluteTime[];
    duration = endTime - startTime;

    (* Remove active job file when done *)
    If[showIntegrands,
        Quiet[DeleteFile[jobFile]];
    ];

    (* Always log completed job; needed for UI progress tracking *)
    logFile = FileNameJoin[{$STJobTrackingDir, "completed_kernel_" <> ToString[$KernelID] <> "_ct_" <> ToString[ctId] <> ".m"}];
    Put[<|
        "Integrand" -> ctIntegrand,
        "Order" -> LRorder,
        "Face" -> faceDirectory,
        "CtId" -> ctId,
        "Duration" -> duration,
        "Success" -> successQ,
        "StartTime" -> startTime,
        "EndTime" -> endTime
    |>, logFile];

    Print["Done"];

    If[Not[DirectoryQ[faceDirectory <> "/partial_results"]],
        CreateDirectory[faceDirectory <> "/partial_results"]
    ];

    (* Write the result file.  New format: an association carrying both
       the integration result and the per-subkernel letter table so the
       main-kernel aggregator can union them.  On failure we still write
       a bare $Failed sentinel (no letter table to carry); the aggregator
       detects both cases. *)
    Put[If[successQ,
           <|"result" -> result, "letterTable" -> localLetterTable|>,
           $Failed],
        faceDirectory <> "/partial_results/result_ct_" <> ToString[ctId] <> ".m"];

    (*  optionally clear HyperIntica memoization caches immediately after each
       integration. When ClearCachesPerIntegrand is True, flush every time; otherwise
       flush only if MemoryInUse[] exceeds ClearCachesMemoryThreshold (given in GB).
       Default Infinity disables threshold-based flushing so cache reuse across
       counter-terms on the same kernel is preserved. *)
    If[OptionValue["ClearCachesPerIntegrand"],
        ForgetAllMemo[]; ClearSystemCache[]
        ,
        If[NumericQ[memThreshold] && MemoryInUse[] > memThreshold * 1024^3,
            ForgetAllMemo[]; ClearSystemCache[]
        ]
    ];

    {ctId, successQ}
]


Options[STLaunchHyperInticaAll] = {
    "NumericsQ" -> False,
    "ShowIntegrands" -> False,
    "ClearCachesPerIntegrand" -> False,
    "ClearCachesMemoryThreshold" -> Infinity,
    "SelectFaces" -> All,
    "ReuseExistingResults" -> True,
    "UIComms" -> None,
    "Integrator" -> "HyperIntica"
};

STLaunchHyperInticaAll[id_:"NP", OptionsPattern[]] := Module[
{NPfilenames, faceFilesToIntegrate, polysAndPairs, xvars, linearlyReducibleOrders, ctOrderDirectory,
integrandsData, ctIntegrands, ctId, order, howManyCt,
ctIdsOutcomes, showIntegrands, selectFaces, reuseQ,
successQ, failed = {}, uiComms},

    If[Not[$KernelSetupQ], STSetupKernel[If[IntegerQ[$STRequestedKernelCount], $STRequestedKernelCount, 3]]];

    uiComms = OptionValue["UIComms"];
    showIntegrands = OptionValue["ShowIntegrands"];
    selectFaces    = OptionValue["SelectFaces"];
    reuseQ         = OptionValue["ReuseExistingResults"];

    (* Clear any old job files *)
    If[showIntegrands && DirectoryQ[$STJobTrackingDir],
        Quiet[DeleteFile /@ FileNames["*.m", $STJobTrackingDir]];
    ];

    NPfilenames = STListDirectoriesNP2[id] // SortBy[#, -#[[2]] &] &;

    (* Apply SelectFaces filter *)
    NPfilenames = STFilterFaceDirectories[NPfilenames, selectFaces];
    If[selectFaces =!= All,
        Print["[SubTropica] SelectFaces: ", Length[NPfilenames], " face director",
              If[Length[NPfilenames] === 1, "y", "ies"], " selected."]
    ];

    (* Apply ReuseExistingResults filter: skip faces that already have a completed result *)
    faceFilesToIntegrate = If[TrueQ[reuseQ],
        Select[NPfilenames, (
            If[STFaceHasResultQ[#[[3]]],
                Print["[SubTropica] Reusing existing result: ", #[[3]]];
                False,
                True
            ]
        ) &],
        NPfilenames
    ];

    OptionValue["NumericsQ"] // Print;
    howManyCt = <||>;
    integrandsData = Join @@ Table[
        ctIntegrands = Get[faceFile <> "/counter_terms_integrands.m"];
        howManyCt[faceFile] = ctIntegrands // Length;
        order = Get[faceFile <> "/bestOrder.m"];
        Table[{faceFile, i, ctIntegrands[[i]], order}, {i, 1, ctIntegrands // Length}]
        , {faceFile, faceFilesToIntegrate[[;; , 3]]}
    ];

    Print["---------------------------------"];
    Print["Launching HyperIntica on " <> ToString[integrandsData // Length] <> " integrands"];

    stUpdateUIComms[uiComms, "JobsTotal", Length[integrandsData]];
    stUpdateUIComms[uiComms, "JobsCompleted", 0];

    DistributeDefinitions[showIntegrands];

    (* Tag each integrand with a 0-based index unique across this ParallelMap.
       Each subkernel uses its JobIndex*Stride as the base offset for its
       Wm[i]/Wp[i] letter allocation, so tables returned from different
       subkernels never collide on the main kernel. When FindRoots is off,
       the Stride is irrelevant (LinearFactors never increments the counter),
       but passing a valid index is still cheap. *)
    Module[{jobStride = $STFindRootsJobStride, indexedData},
        indexedData = MapIndexed[{#1, First[#2] - 1} &, integrandsData];
        (* Capture main-kernel HyperIntica state once and inline its value
           into the ParallelMap closure via With, so each subkernel rebinds
           its flags to these values before calling the integrator.  Using
           a Module-local variable in the closure doesn't round-trip through
           ParallelMap serialization (subkernels see the unbound symbol);
           With forces the Association to be substituted literally. *)
        With[{
                hyperSnap = stBuildHyperSnapshot[],
                numericsQ = OptionValue["NumericsQ"],
                showInt   = showIntegrands,
                clearCpi  = OptionValue["ClearCachesPerIntegrand"],
                clearMT   = OptionValue["ClearCachesMemoryThreshold"],
                integrat  = OptionValue["Integrator"],
                stride    = jobStride
            },
            ctIdsOutcomes = ParallelMap[
                STLaunchHyperInticaAllKernelIntegrator[#[[1]],
                    "NumericsQ" -> numericsQ,
                    "ShowIntegrands" -> showInt,
                    "ClearCachesPerIntegrand" -> clearCpi,
                    "ClearCachesMemoryThreshold" -> clearMT,
                    "Integrator" -> integrat,
                    "JobIndex"      -> #[[2]],
                    "JobStride"     -> stride,
                    "HyperSnapshot" -> hyperSnap
                ] &,
                indexedData,
                Method -> "FinestGrained"
            ];
        ];
    ];

    Print[" integrations completed, formatting results..."];

    (* Only assemble result.m for faces that were actually integrated this run;
       reused faces already have their result.m from a previous run.  If any
       per-counter-term file contains $Failed (the sentinel written when the
       integrator rejected an algebraic-root integrand it couldn't handle),
       propagate $Failed as the face's result instead of wrapping it in
       ZeroInfPeriod[$Failed] via the rewrite rule.

       Per-CT files now come in two formats:
         - legacy:  a bare expression (pre-parallel-safe writer) or $Failed
         - new:     <|"result" -> expr, "letterTable" -> <|...|>|> or $Failed
       We extract the result scalar from both; in the new format we also
       union the per-subkernel letter table into the main-kernel
       $HyperAlgebraicLetterTable so SimplifyWithVieta and
       CanonicalizeAlgebraicLetters see a complete map.  Letter indices
       allocated under different JobIndex values are guaranteed disjoint
       by the stride, so Join-union is safe. *)
    Do[
        Module[{rawCt, ctResults, ctLetterTable, failedQ},
          rawCt = Table[Get[faceFile <> "/partial_results/result_ct_" <> ToString[i] <> ".m"],
                        {i, howManyCt[faceFile]}];
          ctResults = Map[
              Which[
                  # === $Failed, $Failed,
                  Head[#] === Association && KeyExistsQ[#, "result"], #["result"],
                  True, #
              ] &, rawCt];
          ctLetterTable = Join @@ Map[
              If[Head[#] === Association && KeyExistsQ[#, "letterTable"],
                 #["letterTable"], <||>] &, rawCt];
          If[Length[ctLetterTable] > 0,
              HyperIntica`$HyperAlgebraicLetterTable = Join[
                  HyperIntica`$HyperAlgebraicLetterTable, ctLetterTable];
              If[IntegerQ[HyperIntica`$HyperAlgebraicLetterCounter],
                  HyperIntica`$HyperAlgebraicLetterCounter = Max[
                      HyperIntica`$HyperAlgebraicLetterCounter,
                      Max[Keys[ctLetterTable]]]]];
          failedQ = MemberQ[ctResults, $Failed] || AnyTrue[ctResults, (!FreeQ[#, $Failed]) &];
          Put[If[failedQ, $Failed,
                  ctResults //. ZeroInfPeriod :> ZeroInfPeriodAsMpl //. mzv -> zeta],
              faceFile <> "/result.m"];
          Put[!failedQ, faceFile <> "/successQ.m"]
        ]
        , {faceFile, faceFilesToIntegrate[[;; , 3]]}
    ];

    Print[" All done! "];
    Print[" Problems with " <> ToString[Count[ctIdsOutcomes, $Failed]] <> " face integrands!"];
]


(* ::Subsection::Closed:: *)
(*Find LR Subtraction Gauges*)


(* ::Text:: *)
(*Given an integrand with the geometric property and divergent fan \Sigma^{div}, find a superset of the LR gauges. (\[LeftArrow] cannot be done?)*)


(* ::Text:: *)
(*Given an integrand with the geometric property, a choice of v-functions, and the divergent fan \Sigma^{div}, find LR gauges.*)


(* ::Text:: *)
(*Algorithm: *)
(*Step 1: List all possible gauges. *)
(*Step 2: Compute Expansion*)
(*Step 3: Begin Propagation of LR orders*)


(* ::Text:: *)
(*Propagation of LR orders:*)
(*Pre-process: Projectivize/Un-gauge I_\emptyset, call it I_\emptyset again for simplicity.*)
(*Step 1: For each counter-term: enumerate all LR orders (all up to 1 variable). *)
(*Step 2: Intersect the common LR orders *)


(* ::Text:: *)
(*Propagation of LR orders:*)
(*Pre-process: Projectivize/Un-gauge I_\emptyset, call it I_\emptyset again for simplicity.*)
(*Step 1: Enumerate all LR orders (all up to 1 variable). *)


STIIFindGauges[integrand_,xvars_,coeffs_,trData_,divFaces_]:=Module[
{Js},
	(*J space*)
	Js=Table[
		fc-> Cases[Subsets[xvars//Length//Range,{fc//Length}],I_/;Not[Det[trData["rays"][[fc,I]]]===0]:> xvars[[ Complement[Range[xvars//Length],I] ]] ]
		,
		{fc,divFaces[[2;;]]}
	];
	Js=Join[{divFaces[[1]]-> Range[xvars//Length]},Js];
	Js=Association@@Js
]


(* ::Section::Closed:: *)
(*End*)


(*End[];*)
EndPackage[];
(*
On[General::shdw];
(* EndPackage[] triggers a SetFunctionInformation FE packet that the front
   end cannot process (~3700 public symbols). Suppress FE communication
   just for the EndPackage[] call by temporarily nulling $FrontEnd. *)
Block[{$FrontEnd = Null},
  Quiet[EndPackage[], {FrontEndObject::notavail}]];*)
