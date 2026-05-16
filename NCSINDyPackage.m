(* ::Package:: *)

BeginPackage["NCSINDyPackage`"]
(*% Copyright 2026, All Rights Reserved
% Code by Hao Meng
% For Paper, "A noise-corrected data-driven approach for 
nonlinear dynamical system modeling in frequency domain"
% by Jian Tang, Jiahui Xu, Hao Meng, 
Zhihong Zhang, Huidong Xu, Xiaoyan Xiong, Zhihua Wang*)


SVATransform::usage="Convert between Fourier series coefficients of displacement, velocity, and acceleration"

PolynomialTerms::usage="used to Generate all PowerOrder-rank basis functions in the model library"

ModelLibFun::usage="used to generate a power series model library with orders from 1 to PowerOrder"

NonLSLM::usage="Levenberg-Marquardt solver for nonlinear least-squares minimization."

NonLinearSparsifyDynamics::usage="LM-STLS algorithm for nonlinear sparse regression."

AddNoiseFun::usage="Add Gaussian white noise with a specified signal-to-noise ratio (in dB) to a given signal."


Begin["`Private`"]
(* 
Function: SVATransform
Purpose: Convert between Fourier series coefficients of displacement, velocity, and acceleration
Input:
  - inputCoeffs: Input Fourier series coefficients (including DC component), ordered from low to high frequency
    (structure: {DC, cos1, sin1, cos2, sin2, ..., cosN, sinN})
  - freqUnit: Frequency unit ("Hz" or "rad/s")
  - baseFreq: Fundamental frequency (\[Omega]) used in Fourier series expansion
  - operation: Transformation type, where:
    "s2v" = displacement \[RightArrow] velocity
    "s2a" = displacement \[RightArrow] acceleration
    "v2a" = velocity \[RightArrow] acceleration
    "v2s" = velocity \[RightArrow] displacement
    "a2s" = acceleration \[RightArrow] displacement
    "a2v" = acceleration \[RightArrow] velocity
Output:
  - Output Fourier series coefficients in the same structure as inputCoeffs
*)
SVATransform[inputCoeffs_, freqUnit_, baseFreq_, operation_] := 
  Module[{cosCoeffs, sinCoeffs, maxOrder, freqScalars, outputCos, outputSin},
    (* Calculate maximum harmonic order from input length *)
    maxOrder = (Length[inputCoeffs] - 1) / 2;
    
    (* Extract cosine and sine coefficients (excluding DC component) *)
    cosCoeffs = inputCoeffs[[2 ;; -1]][[;;;; 2]];    (* Cosine terms: cos1, cos2, ..., cosN *)
    sinCoeffs = inputCoeffs[[2 ;; -1]][[2;;;; 2]];   (* Sine terms: sin1, sin2, ..., sinN *)
    
    (* Generate frequency scalars (n\[Omega] or n2\[Pi]f) based on unit *)
    freqScalars = Range[maxOrder];  (* Harmonic numbers: 1, 2, ..., N *)
    freqScalars = Which[
      freqUnit === "Hz", 2 \[Pi] baseFreq freqScalars,  (* Convert Hz to angular frequency: n\[CenterDot]2\[Pi]f\:2080 *)
      freqUnit === "rad/s", baseFreq freqScalars,    (* Use rad/s directly: n\[CenterDot]\[Omega]\:2080 *)
      True, MessageDialog["Invalid frequency unit. Use \"Hz\" or \"rad/s\""]; 
             Return[$Failed]
    ];
    
    (* Perform coefficient transformation based on operation type *)
    {outputCos, outputSin} = Switch[
      operation,
      "s2v", {freqScalars * sinCoeffs, -freqScalars * cosCoeffs},    (* Displacement to velocity *)
      "s2a", {-freqScalars^2 * cosCoeffs, -freqScalars^2 * sinCoeffs},  (* Displacement to acceleration *)
      "v2a", {freqScalars * sinCoeffs, -freqScalars * cosCoeffs},    (* Velocity to acceleration *)
      "v2s", {-sinCoeffs / freqScalars, cosCoeffs / freqScalars},    (* Velocity to displacement *)
      "a2s", {-cosCoeffs / freqScalars^2, -sinCoeffs / freqScalars^2},  (* Acceleration to displacement *)
      "a2v", {-sinCoeffs / freqScalars, cosCoeffs / freqScalars},    (* Acceleration to velocity *)
      _, MessageDialog["Invalid operation. Use: \"s2v\";\"s2a\";\"v2a\";\"v2s\";\"a2s\";\"a2v\""]; 
           Return[$Failed]
    ];
    
    (* Format output: {DC, cos1, sin1, cos2, sin2, ...} with DC component = 0 *)
    Join[{0}, Flatten[Table[{outputCos[[j]], outputSin[[j]]}, {j, Length[outputCos]}]]]
  ]


(* 
Purpose:  used to Generate all PowerOrder-rank basis functions in the model library
Input:
  - Vars: List of variables (e.g., {x, x'} are the displacement and velocity used to generate basis
   functions)
  - PowerOrder: Integer specifying the degree of model to generate (e.g., 2 for quadratic terms)
  - FourierOrder: (Optional) Maximum harmonic order for Fourier series conversion (None = no conversion)
Output:
  - If FourierOrder = None: List of pure Terms (e.g., x\.b2, x*x', x'*x' for PowerOrder=2)
  - If FourierOrder \[NotEqual] None: The related Expressions in the frequency domain, each with 
    structure {constant_term, cos(\[Omega]t)_coeff, sin(\[Omega]t)_coeff, ..., cos(n\[Omega]t)_coeff, sin(n\[Omega]t)_coeff}
*)
PolynomialTerms[Vars_,PowerOrder_,FourierOrder_:None]:=
  Module[{numVars,tempVars,polyExpr,expandedPoly,monomialList,coeffs,pureTerms,temp01,j1,temp02,outPut},
  numVars=Length[Vars]; (* Number of input variables *)
  
  (* Define temporary variables for polynomial expansion (x1, x2, ..., xn) *)
  tempVars=If[FourierOrder=!=None,Array[x,numVars],Vars];
  
  (* Generate (x1 + x2 + ... + xn)^PowerOrder and expand into monomials *)
  polyExpr=Total[tempVars]^PowerOrder;   (* Sum of variables raised to PowerOrder *)
  expandedPoly=Expand[polyExpr];         (* Expand into individual monomials *)
  
  (* Extract monomials (handle single-variable case specially to avoid list issues) *)
  monomialList=If[numVars==1,{expandedPoly},List@@expandedPoly];
  
  (* Extract coefficients of each monomial (set all variables to 1 to isolate coefficients) *)
  coeffs=monomialList/.Thread[tempVars->1];
  
  (* Remove coefficients to retain pure monomial structures (e.g., 2x1x2 \[RightArrow] x1x2) *)
  pureTerms=monomialList/coeffs;
  
  (* If FourierOrder is specified, convert pureTerms to their frequency domain Expressions *)
  If[FourierOrder=!=None,
    temp01=TrigReduce[pureTerms/.Thread[tempVars->Flatten[Vars.{Join[{1},Flatten[Array[{Cos[# \[CapitalOmega] t],
            Sin[# \[CapitalOmega] t]}&,(Dimensions[Vars][[2]]-1)/2]]]}\[Transpose]]]];
      
    For[j1=1,j1<=Length[temp01],j1++,
      temp02[j1]=Join[{temp01[[j1]]}/.{Sin[a_.*t]:>0,Cos[a_.*t]:>0},Coefficient[temp01[[j1]],
                  Flatten[Array[{Cos[# \[CapitalOmega] t],Sin[# \[CapitalOmega] t]}&,FourierOrder]]]]
    ];
    outPut=Array[temp02[#]&,Length[temp01]]\[Transpose],
    outPut=pureTerms
  ];
outPut]


(* 
Purpose: used to generate a power series model library with orders from 1 to PowerOrder
Input:
  - Vars: List of variables (same as in PolynomialTerms, e.g., {x, x'})
  - MaxPowerOrder: Maximum degree of model to include (e.g., 2 includes 1st and 2nd order terms)
  - FourierOrder: (Optional) Maximum harmonic order for Fourier series conversion (None = no conversion)
Output:
  - Model library of basis functions with orders from 1 to MaxPowerOrder:
    - If FourierOrder = None: Terms included in Model library (e.g., {x, x', x\.b2, x*x', x'*x'} for 
    MaxPowerOrder=2)
    - If FourierOrder \[NotEqual] None: Model library in frequency domain
*)
ModelLibFun[Vars_,MaxPowerOrder_,FourierOrder_:None]:=
  Module[{PowerOrder,Temp01,OutputModelLib},
  OutputModelLib=If[FourierOrder=!=None,{{}},{}];
  For[PowerOrder=1,PowerOrder<=MaxPowerOrder,PowerOrder++,
    Temp01[PowerOrder]=PolynomialTerms[Vars,PowerOrder,FourierOrder];
    OutputModelLib=If[FourierOrder=!=None,
      Join[OutputModelLib,Temp01[PowerOrder],2],Join[OutputModelLib,Temp01[PowerOrder]]];
  ];
OutputModelLib]


(* 
  Levenberg-Marquardt solver for nonlinear least-squares minimization.
  Minimizes the objective function F[\[Beta]] = 1/2 ||r(\[Beta])||\.b2, where r(\[Beta]) is the residual vector.
  
  Inputs:
    - jacobianFun: Function handle for Jacobian matrix J[\[Beta]] (m\[Times]n, m \[GreaterEqual] n)
    - residualFun: Function handle for residual vector r[\[Beta]] (length m)
    - variables: List of parameter symbols {\[Beta]\:2081, \[Beta]\:2082, ..., Subscript[\[Beta], n]}
    - initGuess: Initial parameter vector \[Beta]\:2080 (length n)
    - convTol: Convergence tolerance (stopping criterion for parameter update norm)
    - maxIter: Maximum number of iterations
    - initDamping: Initial damping parameter \[Mu]\:2080 (default: 0.01)
    - dampingFactor: Initial factor for increasing damping (default: 10)
  
  Outputs:  {solCurr, convFlag, iter, stepNorm, \[Mu], solHistory}
    - solCurr: Optimal parameter vector \[Beta]*
    - convFlag: Convergence flag (1 = converged to tolerance; 0 = max iterations reached)
    - iterCount: Number of iterations executed
    - stepNorm: Maximum norm of the final parameter update (for convergence validation)
    - \[Mu]: Final value of the damping parameter
    - solHistory: History of variable estimates across iterations
*)
NonLSLM[jacobianFun_, residualFun_, variables_, initGuess_, convTol_, maxIter_, initDamping_: 0.01, dampingFactor_: 10] := 
 Module[{solHistory={}, jFun = jacobianFun, rFun = residualFun, var = variables, solCurr ,
          \[Mu] = initDamping, \[Nu] = dampingFactor, iterCount = 0, 
         AA, RR, solStep, convFlag, newCost, curCost, \[Rho], stepAccepted, stepNorm},

  (* Step 1: Initialize parameters and initial residual cost *)
  solCurr = N[initGuess];
  curCost = Total[(rFun /. Thread[var -> solCurr])^2];
  
  (* Step 2: Iterative LM optimization loop *)
  While[iterCount < maxIter,
    iterCount++;
    solHistory=AppendTo[solHistory,solCurr];  (* Store current parameter estimate *)
    AA = N[jFun /. Thread[var -> solCurr]];  (* evaluate Jacobian matrix *)
    RR = N[-(rFun /. Thread[var -> solCurr])]; (* residual cost vector *)
    
    (* Compute LM update step *)
    solStep = LeastSquares[Transpose[AA].AA + \[Mu] IdentityMatrix[Length[solCurr]], Transpose[AA].RR];
    
    (* Evaluate cost at candidate parameters (solCurr + solStep) *)
    newCost = Total[(rFun /. Thread[var -> (solCurr + solStep)])^2];
    
    (* Calculate reduction ratio \[Rho] to assess step quality *)
    \[Rho] = (curCost - newCost) / (solStep.(\[Mu] solStep + Transpose[AA].RR));

    (* Adaptively adjust damping and update parameters *)
    If[\[Rho] > 0,
      stepAccepted = True;
      solCurr=solCurr + solStep;
      curCost = newCost;
      \[Mu] *= Max[1/3, 1 - (2\[Rho] - 1)^3];
      \[Nu] = 10,
      (* Else *)
      stepAccepted = False;
      \[Mu] *= \[Nu];
      \[Nu] *= 2
    ];
    
    (* Check convergence (small parameter update) *)
    stepNorm=Max[Norm[solStep]];
    If[stepNorm < convTol, Break[]];
  ];
  
    (* Step 3: Set convergence flag and prepare outputs *)
    convFlag = If[iterCount >= maxIter, 0, 1];
  
    (* Return results *)
    {solCurr, convFlag, iterCount, stepNorm, \[Mu], solHistory}
]


(* 
  LM-STLS algorithm for nonlinear sparse regression.
  Alternates between STLS thresholding (sparsity promoting) and LM-based refinement to minimize 
  ||r(\[Beta])||\.b2 while preserving critical regressors.
  
  Inputs:
    - jacobianFun: Function handle for Jacobian matrix J[\[Beta]]
    - residualFun: Function handle for residual vector r[\[Beta]]
    - variables: List of variable symbols {\[Beta]\:2081, \[Beta]\:2082, ..., Subscript[\[Beta], p]}
    - initGuess: Initial variable vector \[Beta]\:2080
    - convTol: Convergence tolerance for LM subproblems (passed to NonLSLM)
    - maxLMIter: Max iterations for LM solver per refinement (passed to NonLSLM)
    - lambda: Threshold lambda for sparsity (parameters with |\[Beta]\:1d62| < lambda are zeroed)
    - maxSparseIter: Maximum number of sparsification cycles (default: 10)
    - initbigIndices: Indices of parameters to retain regardless of threshold (physical constraints)
  
  Outputs:
    - paramHistory: Full history of parameter estimates across all cycles/iterations
    - iterCountHistory: Number of LM iterations used in each sparsification cycle
*)
NonLinearSparsifyDynamics[jacobianFun_,residualFun_,variables_,initGuess_,convTol_,
                     maxLMIter_,lambda_,maxSparseIter_: 10,initbigIndices_]:=
  Module[{tempLM,paramHistory={},iterCountHistory={},solHistory,smallIndices,largeIndices,
  substitutionRule,reducedJacobian,reducedResidual,reducedVars,reducedInitGuess,tempsolHistory,
ZeroMatrix,k1},
   (* Step 1: Initial dense solution via LM algorithm *)
  tempLM=NonLSLM[jacobianFun,residualFun,variables,initGuess,convTol,maxLMIter];
  solHistory[0]=tempLM[[1]]; (* Store initial parameter estimate (solution) *)
  
   (* Initialize history with results from the initial LM run *)
  paramHistory=Join[paramHistory,tempLM[[-1]]]; (* Append initial LM iteration history *)
  iterCountHistory=AppendTo[iterCountHistory,tempLM[[3]]]; (* Track LM iterations used *)
  
  (* Step 2: Iterative sparsification cycles *)
  For[k1=1,k1<=maxSparseIter,k1++,
    solHistory[k1-1]=solHistory[k1-1]//Flatten;
    
    (* Identify small parameters to zero (excluding initRetainIndices) using STLS *)
    smallIndices=Position[solHistory[k1-1],x_/;Abs[x]<lambda]//Flatten;
    smallIndices=Complement[smallIndices,initbigIndices];
    
    (* Apply thresholding: zero out small parameters *)
    solHistory[k1]=ReplacePart[solHistory[k1-1],Thread[smallIndices->0]];
    solHistory[k1]=If[Length[Dimensions[solHistory[k1]]]!=1,solHistory[k1],{solHistory[k1]}\[Transpose]];
    
    (* Identify retained parameters (larger than threshold value lambda) *)
    largeIndices=Union[Position[solHistory[k1][[All,1]],x_/;Abs[x]>=lambda]//Flatten,initbigIndices];
    
    (* Reduce problem to retained variables *)
    substitutionRule=Thread[variables[[smallIndices//Flatten]]->0];
    reducedJacobian=jacobianFun[[All,Flatten[largeIndices]]]/.substitutionRule;
    reducedResidual=If[Length[Dimensions[residualFun]]==1,residualFun/.substitutionRule,
                      residualFun[[All,1]]/.substitutionRule];
    reducedVars=variables[[largeIndices//Flatten]];
    reducedInitGuess=solHistory[k1][[largeIndices//Flatten]]//Flatten;
    
    (* Refine retained variables using LM *)
    tempLM=NonLSLM[reducedJacobian,reducedResidual,reducedVars,reducedInitGuess,convTol,maxLMIter];
    
    (* Update variables history (pad excluded parameters with zeros) *)
    tempsolHistory=tempLM[[-1]];ZeroMatrix=ConstantArray[0,{Dimensions[tempsolHistory][[1]],Length[variables]}];
    Table[ZeroMatrix[[All,largeIndices[[j1]]]]=tempsolHistory[[All,j1]],{j1,Length[largeIndices]}];
    paramHistory=Join[paramHistory,ZeroMatrix];
    
    (* Track LM iterations and update current solution *)
    iterCountHistory=AppendTo[iterCountHistory,tempLM[[3]]];
    solHistory[k1]=ReplacePart[solHistory[k1],Thread[Thread[{Flatten[largeIndices],1}]->tempLM[[1]]]];
    Print[{k1,Count[Flatten[solHistory[k1]],x_/;x!=0]}];];
    
    (* Return full history and iteration counts *)
    {paramHistory,iterCountHistory}
]


(* 
  Add Gaussian white noise with specified SNR (in dB) to a given signal.
  Generates noise with desired power based on input signal power and target SNR.
  
  Inputs:
    - Sig: Original clean signal vector (1D array)
    - MeanNoise: Mean value of the Gaussian white noise
    - SNRdb: Target signal-to-noise ratio in decibels (dB)
  
  Outputs:
    - noise: Generated Gaussian white noise vector with specified statistical properties
*)
AddNoiseFun[Sig_,MeanNoise_,SNRdb_]:=Module[{n,PSignal,SNR,Pnoise,SigmaNoise,noise},
  (* Step 1: Obtain signal length *)
n=Length[Sig];
(* Step 2: Compute average power of the input signal *)
PSignal=Norm[Sig]^2/Length[Sig];
 (* Step 3: Convert SNR from decibels (dB) to linear scale *)
SNR=10^(SNRdb/10);
 (* Step 4: Calculate noise power based on signal power and SNR *)
Pnoise=PSignal/SNR;
 (* Step 5: Compute noise standard deviation from noise power *)
SigmaNoise=Sqrt[Pnoise];
(* Step 6: Generate Gaussian white noise with specified mean and standard deviation *)
noise=SigmaNoise RandomReal[NormalDistribution[MeanNoise,1],n];
(* Return generated noise *)
noise]


End[]

EndPackage[]
