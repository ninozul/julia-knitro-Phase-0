# model for ARPA-E Competition Phase 0
#
# getting constraint violations for real and reactive power balance at all buses
# in the base case and all contingency cases.
#
# one possible problem is in the cos() and sin() expressions.
# I believe these expressions do not match the documentation exactly
# the documentation has the argument as
#
#   theta_to - theta_from + theta_transformer
#
# together with the +/- signs on the cos and sin terms that this model uses,
# I believe this model would be correct if the term
#
#   theta_transformer
#
# carried a factor of -1, i.e. a (-) sign rather than a (+) sign.
# not sure though as this depends on the data input process.
# Not correcting this yet.
# Better to see what the evaluation code comes up with after fixing
# the following error.
# it is possible that the code is correct, depending on the data input process.
#
# another possible problem is in the output of the voltage angles.
# it appears that the maximum angle difference in this model comes out as
# about 3 degrees (in solution2.txt), while the maximum angle difference in the
# gams model is 10 or 15 degrees.
# SOLVED: the issue here is the output to solution2.txt uses v and v0,
# rather than theta and theta0. This has been corrected.
#
# these observations are all without any pv/pq switching in the model.

using Complementarity;
#using Ipopt;
using KNITRO;
#using Ipopt, Complementarity;

# Global: contDList

function buildMod(fData,uData, contDList)
  # Input:
  # fData - grid data with all parameters
  # uData - contingency data

  # readin the data from the structure
  baseMVA = fData.baseMVA;
  bList = fData.busList;
  bData = fData.busDList;

  gList = fData.genList;
  gData = fData.genDList;

  brList = fData.brList;
  brData = fData.brDList;

  S = uData.contList;
  contData = uData.contDList;

  # set up the model
  #mp = Model(solver = IpoptSolver(print_level=0));
  mp = Model(solver = KnitroSolver(KTR_PARAM_OUTLEV=2,  # default is 2
   			feastol=2.25e-9, 
   			#feastol_abs=1e-2,
   			opttol=1e-4, 
   			cg_maxit=10,   # formerly maxcgit
   			maxit=400,
   			ftol=1e-4, 
   			#ftol_iters=3, 
   			#pivot=1e-12,
        #ms_enable=1,
        #ms_maxsolves=5,
        #ms_maxbndrange=10,
        #ms_terminate=1,
        #par_numthreads=24,
        #par_concurrent_evals=0,
   			maxtime_real=3600)); 
  
  # create the variables for the base case
  @variable(mp,bData[i].Vmin <= v0[i in bList] <= bData[i].Vmax);
  @variable(mp,gData[l].Pmin <= sp0[l in gList] <= gData[l].Pmax);
  @variable(mp,gData[l].Qmin <= sq0[l in gList] <= gData[l].Qmax);
  @variable(mp,p0[k in brList]);
  @variable(mp,q0[k in brList]);
  @variable(mp,psh0[i in bList]);
  @variable(mp,qsh0[i in bList]);
  @variable(mp,theta0[i in bList]);

  # create the constraints for the base case
  @NLconstraint(mp,flowBound0[k in brList],p0[k]^2 + q0[k]^2 <= brData[k].t^2);
  @NLconstraint(mp,pShunt0[i in bList],psh0[i] == bData[i].gsh*v0[i]^2);
  @NLconstraint(mp,qShunt0[i in bList],qsh0[i] == -bData[i].bsh*v0[i]^2);
  @NLconstraint(mp,pFlow0[k in brList;brData[k].zeroImpe == false], p0[k] == brData[k].g/(brData[k].tauprime^2)*v0[brData[k].From]^2
               - brData[k].g/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*cos(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               - brData[k].b/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*sin(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               );
  @NLconstraint(mp,qFlow0[k in brList;brData[k].zeroImpe == false], q0[k] == (-brData[k].b - brData[k].bc/2)/(brData[k].tauprime^2)*v0[brData[k].From]^2
               + brData[k].b/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*cos(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               - brData[k].g/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*sin(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               );
  @constraint(mp,zeroP0[k in brList;brData[k].zeroImpe == true], p0[k] == -p0[brData[k].revID]);
  @NLconstraint(mp,zeroQ0[k in brList;brData[k].zeroImpe == true], q0[k] + q0[brData[k].revID] == -brData[k].bc*v0[brData[k].To]^2);
  @constraint(mp,zeroV0[k in brList;brData[k].zeroImpe == true], v0[brData[k].To] == v0[brData[k].From]/brData[k].tau);
  @constraint(mp,zerotheta0[k in brList;brData[k].zeroImpe == true], theta0[brData[k].To] == theta0[brData[k].From] - brData[k].thetatr);
  @constraint(mp,pBalance0[i in bList], sum(sp0[l] for l in bData[i].gen) == psh0[i] + bData[i].Pd + sum(p0[k] for k in brList if (brData[k].From == i)));
  @constraint(mp,qBalance0[i in bList], sum(sq0[l] for l in bData[i].gen) == qsh0[i] + bData[i].Qd + sum(q0[k] for k in brList if (brData[k].From == i)));

  # create the variables for the contingency cases
  @variable(mp,bData[i].Vmin <= v[i in bList, s in S] <= bData[i].Vmax);
  @variable(mp,gData[l].Pmin <= sp[l in gList, s in S;(!(l in contData[s].Loc))] <= gData[l].Pmax);
  @variable(mp,gData[l].Qmin <= sq[l in gList, s in S;(!(l in contData[s].Loc))] <= gData[l].Qmax);
  @variable(mp,p[k in brList, s in S; !(k in contDList[s].Loc)]);
  @variable(mp,q[k in brList, s in S; !(k in contDList[s].Loc)]);
  @variable(mp,psh[i in bList, s in S]);
  @variable(mp,qsh[i in bList, s in S]);
  @variable(mp,theta[i in bList, s in S]);
  @variable(mp,pdelta[s in S]);

  # create the constraints for the contingency cases
  @NLconstraint(mp,flowBoundS[k in brList, s in S; !(k in contDList[s].Loc)],p[k,s]^2 + q[k,s]^2 <= brData[k].t^2);
  @NLconstraint(mp,pShuntS[i in bList, s in S],psh[i,s] == bData[i].gsh*v[i,s]^2);
  @NLconstraint(mp,qShuntS[i in bList, s in S],qsh[i,s] == -bData[i].bsh*v[i,s]^2);
  @NLconstraint(mp,pFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], p[k,s] == brData[k].g/(brData[k].tauprime^2)*v[brData[k].From,s]^2
               - brData[k].g/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*cos(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               - brData[k].b/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*sin(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               );
  @NLconstraint(mp,qFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], q[k,s] == (-brData[k].b - brData[k].bc/2)/(brData[k].tauprime^2)*v[brData[k].From,s]^2
               + brData[k].b/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*cos(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               - brData[k].g/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*sin(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               );
  @constraint(mp,zeroPS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], p[k,s] == -p[brData[k].revID,s]);
  @NLconstraint(mp,zeroQS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], q[k,s] + q[brData[k].revID,s] == -brData[k].bc*v[brData[k].To,s]^2);
  @constraint(mp,zeroVS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], v[brData[k].To,s] == v[brData[k].From,s]/brData[k].tau);
  @constraint(mp,zerothetaS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], theta[brData[k].To,s] == theta[brData[k].From,s] - brData[k].thetatr);
  @constraint(mp,pBalanceS[i in bList, s in S], sum(sp[l,s] for l in bData[i].gen if (!(l in contData[s].Loc))) == psh[i,s] + bData[i].Pd + sum(p[k,s] for k in brList if ((brData[k].From == i)&(!(k in contDList[s].Loc)))));
  @constraint(mp,qBalanceS[i in bList, s in S], sum(sq[l,s] for l in bData[i].gen if (!(l in contData[s].Loc))) == qsh[i,s] + bData[i].Qd + sum(q[k,s] for k in brList if ((brData[k].From == i)&(!(k in contDList[s].Loc)))));
  @constraint(mp,recourse[l in gList, s in S; (!(l in contData[s].Loc))], sp[l,s] == sp0[l] + gData[l].alpha*pdelta[s]);

  # create the Complementarity constraint
  @variable(mp, delplus[i in bList, s in S] >= 0);
  @variable(mp, delminus[i in bList, s in S] >= 0);
  @variable(mp, sqplus[i in bList,s in S] >= 0);
  @variable(mp, sqminus[i in bList,s in S] >= 0);
  @constraint(mp,vConstr1[i in bList,s in S], v0[i] - v[i,s] <= delplus[i,s]);
  @constraint(mp,vConstr2[i in bList,s in S], v[i,s] - v0[i] <= delminus[i,s]);
  @constraint(mp, spplusConstr[i in bList,s in S], sqplus[i,s] == sum(gData[l].Qmax - sq[l,s] for l in bData[i].gen));
  @constraint(mp, spminusConstr[i in bList,s in S], sqminus[i,s] == sum(sq[l,s] - gData[l].Qmin for l in bData[i].gen));
  for s in S
    for i in bList
      @complements(mp,0 <= delplus[i,s], sqplus[i,s] >= 0);
      @complements(mp,0 <= delminus[i,s], sqminus[i,s] >= 0);
    end
  end

#  # add the contingency effect as constraints
#  for s in S
#    # for each contingency
#    # first identify the type of the contingency
#    if contData[s].Type == "B"
#      # if it is a branch contingency
#      for k in contData[s].Loc
#        @constraint(mp, p[k,s] == 0.0);       # the real power flow through the branch is 0
#        @constraint(mp, q[k,s] == 0.0);       # the reactive power flow through the branch is 0
#      end
#    elseif contData[s].Type == "T"
#      # if it is a branch contingency
#      for k in contData[s].Loc
#        @constraint(mp, p[k,s] == 0.0);       # the real power flow through the branch is 0
#        @constraint(mp, q[k,s] == 0.0);       # the reactive power flow through the branch is 0
#      end
#    elseif contData[s].Type == "G"
#      # if it is a generator contingency
#      for l in contData[s].Loc
#        @constraint(mp, sp[l,s] == 0.0);
#        @constraint(mp, sq[l,s] == 0.0);
#      end
#    end
#  end

  # build the objective function
  @NLobjective(mp, Min, sum(sum(gData[l].cParams[n]*(sp0[l]*fData.baseMVA)^n for n in gData[l].cn) for l in gList));
  return mp;
end
