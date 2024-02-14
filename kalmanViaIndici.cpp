// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
const double M_LOG2PI = 1.8378770664093454835606594728112352797227949472756;

// [[Rcpp::export]]
List kalman(arma::sp_mat T,
            arma::sp_mat Q,
            arma::sp_mat Z,
            arma::sp_mat H,
            arma::sp_mat Q1,
            arma::sp_mat R,
            const arma::sp_mat a1,
            const arma::sp_mat P1,
            const arma::mat X,
            const arma::mat iZ,
            const arma::mat iH,
            const arma::mat iQ,
            const arma::mat iT,
            const arma::mat iQ1,
            const arma::mat iR,
            const arma::mat Yt,
            const List Smo = "no"
              
) {
  const arma::uword n = Yt.n_cols; // n. of observations
  const arma::uword m = Yt.n_rows; // n. of observed variables
  const arma::uword p = T.n_rows;  // n. of state variables
  
  const arma::uword rH = H.n_rows;// righe matrice H
  const arma::uword rQ1 = Q1.n_rows;// righe matrice Q1
  
  double loglik=0;
  arma::sp_mat att(p, n);
  arma::sp_mat at(p, n+1);
  arma::sp_mat Ythat(m, n);
  arma::sp_mat v(m, n);
  arma::cube Ptt(p, p, n);
  arma::cube Pt(p, p, n+1);
  arma::cube F(m, m, n);
  arma::cube K(p, m, n);
  arma::sp_mat PZt(p, m);

  arma::sp_mat iF(m, m);
  arma::uvec notna;
  arma::colvec v0;
  arma::sp_mat PZ0t;
  at.col(0) = a1;
  Pt.slice(0) = P1;
  
  // smoother
  arma::sp_mat r(p, n);
  arma::sp_mat alphahat(p,n);
  arma::cube N(p, p, n);
  arma::cube Vt(p, p, n);
  arma::cube Lt(p, p, n);
  arma::sp_mat ZiF(p, p);
  arma::cube Kt(p, m, n);
  
  arma::colvec r0;
  arma::sp_mat N0(p,p);
  
  arma::sp_mat VV(m,p);
  arma::uword rz;
  arma::uword cz;

  arma::sword a = 1;
  arma::sword b = 1;
  arma::sword c = 1;
  arma::sword d = 1;
  arma::sword e = 1;
  arma::sword f = 1;
  
  arma::uword dd;
  
  // disurbance
  arma::sp_mat ut(m,n);
  arma::sp_mat epshat(rH,n);
  arma::sp_mat etahat(rQ1,n);
  arma::cube D(m,m,n);
  arma::sp_mat V_eps(m,n);
  arma::cube V_eta(rQ1,rQ1,n);
  
  for (arma::uword t=0; t < n; ++t) {
    dd = 0;
    
    if (!(iZ.n_rows == 2)) {
      while (iZ(a,0) == t){
        rz = iZ(a,2);
        cz = iZ(a,1);
        Z(rz,cz) = X(dd,t);
        dd++;
        a++;
        }
    }
    
    if (!(iT.n_rows == 2)) {
      while (iT(c,0) == t){
        rz = iT(c,2);
        cz = iT(c,1);
        T(rz,cz) = X(dd,t);
        dd++;
        c++;
      }
    }
    
    if (!(iQ.n_rows == 2)) {
      while (iQ(d,0) == t){
        rz = iQ(d,2);
        cz = iQ(d,1);
        Q(rz,cz) = X(dd,t);
        dd++;
        d++;
      }
    }
    
    if (!(iH.n_rows == 2)) {
      while (iH(b,0) == t){
        rz = iH(b,2);
        cz = iH(b,1);
        H(rz,cz) = X(dd,t);
        dd++;
        b++;
      }
    }
    
    if (!(iQ1.n_rows == 2)) {
      while (iQ1(e,0) == t){
        rz = iQ1(e,2);
        cz = iQ1(e,1);
        Q1(rz,cz) = X(dd,t);
        dd++;
        e++;
      }
    }
    
    if (!(iR.n_rows == 2)) {
      while (iR(f,0) == t){
        rz = iR(f,2);
        cz = iR(f,1);
        R(rz,cz) = X(dd,t);
        dd++;
        f++;
      }
    }

    notna = find_finite(Yt.col(t));
    if (notna.n_elem == m) { // no missing obs in y
      // innovation
      Ythat.col(t) = Z*at.col(t);
      v.col(t) = Yt.col(t) - Ythat.col(t);
      PZt = Pt.slice(t)*Z.t();
      F.slice(t) = Z*PZt + H;
      iF = inv_sympd(F.slice(t));
      loglik += (log_det_sympd(F.slice(t)) +
        v.col(t).t()*iF*v.col(t) +
        m*M_LOG2PI).eval().at(0,0);
      // update
      K.slice(t) = PZt * iF;
      att.col(t) = at.col(t) + K.slice(t)*v.col(t);
      Ptt.slice(t) = Pt.slice(t) - PZt*K.slice(t).t();
      // prediction
      at.col(t+1) = T*att.col(t);
      Pt.slice(t+1) = T*Ptt.slice(t)*T.t() + Q;
      
      continue;
    }
    
    if (notna.n_elem == 0) { // all missing obs in y
      // innovation
      Ythat.col(t) = Z*at.col(t);
      v.col(t) = Yt.col(t) - Ythat.col(t);
      PZt = Pt.slice(t)*Z.t();
      F.slice(t) = Z*PZt + H;
      // update
      att.col(t) = at.col(t);
      Ptt.slice(t) = Pt.slice(t);
      // prediction
      at.col(t+1) = T*att.col(t);
      Pt.slice(t+1) = T*Ptt.slice(t)*T.t() + Q;
      
      continue;
    }
    
    // partially missing values in Y
    // innovation
    Ythat.col(t) = Z*at.col(t);
    v.col(t) = Yt.col(t) - Ythat.col(t);
    v0 = arma::colvec(v.col(t))(notna);
    PZt = Pt.slice(t)*Z.t();
    PZ0t = PZt.cols(notna);
    F.slice(t) = Z*PZt + H;
    iF = inv_sympd(F.slice(t)(notna, notna));
    loglik += (log_det_sympd(F.slice(t)(notna, notna)) +
      v0.t()*iF*v0).eval()(0,0) + notna.n_elem*M_LOG2PI;
    // update
    K.slice(t).cols(notna) = PZ0t * iF;
    att.col(t) = at.col(t) + K.slice(t).cols(notna)*v0;
    Ptt.slice(t) = Pt.slice(t) - PZ0t*K.slice(t).cols(notna).t();
    // prediction
    at.col(t+1) = T*att.col(t);
    Pt.slice(t+1) = T*Ptt.slice(t)*T.t() + Q;
    
  }
  
  if (Smo.size() == 1){
    if (Rcpp::as<std::string>(Smo[0]) == "Smoother") {
      a = iZ.n_rows-2;
      b = iT.n_rows-2;
      c = iQ.n_rows-2;
      d = iH.n_rows-2;
      e = iQ1.n_rows-2;
      f = iR.n_rows-2;
      for (arma::uword t=n-1; t>=1; --t) {
        
        dd = X.n_rows - 1;
        
        if (!(iR.n_rows == 2)) {
          while (iR(f,0) == t){
            rz = iR(f,2);
            cz = iR(f,1);
            R(rz,cz) = X(dd,t);
            dd--;
            f--;
          }
        }
        
        if (!(iQ1.n_rows == 2)) {
          while (iQ1(e,0) == t){
            rz = iQ1(e,2);
            cz = iQ1(e,1);
            Q1(rz,cz) = X(dd,t);
            dd--;
            e--;
          }
        }
        
        if (!(iH.n_rows == 2)) {
        while (iH(d,0) == t){
          rz = iH(d,2);
          cz = iH(d,1);
          H(rz,cz) = X(dd,t);
          dd--;
          d--;
        }
       }
        
        if (!(iQ.n_rows == 2)) {
          while (iQ(c,0) == t){
            rz = iQ(c,2);
            cz = iQ(c,1);
            Q(rz,cz) = X(dd,t);
            dd--;
            c--;
          }
        }
        
        if (!(iT.n_rows == 2)) {
          while (iT(b,0) == t){
            rz = iT(b,2);
            cz = iT(b,1);
            T(rz,cz) = X(dd,t);
            dd--;
            b--;
          }
        }
        
        if (!(iZ.n_rows == 2)) {
          while (iZ(a,0) == t){
            rz = iZ(a,2);
            cz = iZ(a,1);
            Z(rz,cz) = X(dd,t);
            dd--;
            a--;
            
          }
        } 
      iF = inv_sympd(F.slice(t));
      Kt.slice(t) = T*Pt.slice(t)*Z.t()*iF;
      Lt.slice(t) = T - Kt.slice(t)*Z;
      ZiF = Z.t()*iF;
      r.col(t-1) = ZiF*v.col(t) + Lt.slice(t).t()*r.col(t); 
      N.slice(t-1) = ZiF*Z + Lt.slice(t).t()*N.slice(t)*Lt.slice(t);
      alphahat.col(t) = at.col(t) + Pt.slice(t)*r.col(t-1);   
      Vt.slice(t) = Pt.slice(t) - Pt.slice(t)*N.slice(t-1)*Pt.slice(t); 
           

    dd = X.n_rows - 1;
    
      if (!(iR.n_rows == 2)) {
        while (iR(f,0) == 0){
          rz = iR(f,2);
          cz = iR(f,1);
          R(rz,cz) = X(dd,0);
          dd--;
          f--;
        }
      }
      
      if (!(iQ1.n_rows == 2)) {
        while (iQ1(e,0) == 0){
          rz = iQ1(e,2);
          cz = iQ1(e,1);
          Q1(rz,cz) = X(dd,0);
          dd--;
          e--;
        }
      }
      
      if (!(iH.n_rows == 2)) {
        while (iH(d,0) == 0){
          rz = iH(d,2);
          cz = iH(d,1);
          H(rz,cz) = X(dd,0);
          dd--;
          d--;
        }
      }
      
      if (!(iQ.n_rows == 2)) {
        while (iQ(c,0) == 0){
          rz = iQ(c,2);
          cz = iQ(c,1);
          Q(rz,cz) = X(dd,0);
          dd--;
          c--;
        }
      }
      
      if (!(iT.n_rows == 2)) {
        while (iT(b,0) == 0){
          rz = iT(b,2);
          cz = iT(b,1);
          T(rz,cz) = X(dd,0);
          dd--;
          b--;
        }
      }
      
      if (!(iZ.n_rows == 2)) {
        while (iZ(a,0) == 0){
          rz = iZ(a,2);
          cz = iZ(a,1);
          Z(rz,cz) = X(dd,0);
          dd--;
          a--;
        }
      } 
    
    iF = inv_sympd(F.slice(0));
    Kt.slice(0) = T*Pt.slice(0)*Z.t()*iF;
    Lt.slice(0) = T - Kt.slice(0)*Z;
    ZiF = Z.t()*iF;
    r0 = ZiF*v.col(0) + Lt.slice(0).t()*r.col(0); 
    N0 = ZiF*Z + Lt.slice(0).t()*N.slice(0)*Lt.slice(0); 
    alphahat.col(0) = at.col(0) + Pt.slice(0)*r0;    
    Vt.slice(0) = Pt.slice(0) - Pt.slice(0)*N0*Pt.slice(0);
  }
    }
    
    if (Rcpp::as<std::string>(Smo[0]) == "Disturbance"){
      a = iZ.n_rows-2;
      b = iT.n_rows-2;
      c = iQ.n_rows-2;
      d = iH.n_rows-2;
      e = iQ1.n_rows-2;
      f = iR.n_rows-2;
      
      for (arma::uword t=n-1; t>=1; --t) {
        
        dd = X.n_rows - 1;
        
        if (!(iR.n_rows == 2)) {
          while (iR(f,0) == t){
            rz = iR(f,2);
            cz = iR(f,1);
            R(rz,cz) = X(dd,t);
            dd--;
            f--;
          }
        }
        
        if (!(iQ1.n_rows == 2)) {
          while (iQ1(e,0) == t){
            rz = iQ1(e,2);
            cz = iQ1(e,1);
            Q1(rz,cz) = X(dd,t);
            dd--;
            e--;
          }
        }
        
        if (!(iH.n_rows == 2)) {
          while (iH(d,0) == t){
            rz = iH(d,2);
            cz = iH(d,1);
            H(rz,cz) = X(dd,t);
            dd--;
            d--;
          }
        }
        
        if (!(iQ.n_rows == 2)) {
          while (iQ(c,0) == t){
            rz = iQ(c,2);
            cz = iQ(c,1);
            Q(rz,cz) = X(dd,t);
            dd--;
            c--;
          }
        }
        
        if (!(iT.n_rows == 2)) {
          while (iT(b,0) == t){
            rz = iT(b,2);
            cz = iT(b,1);
            T(rz,cz) = X(dd,t);
            dd--;
            b--;
          }
        }
        
        if (!(iZ.n_rows == 2)) {
          while (iZ(a,0) == t){
            rz = iZ(a,2);
            cz = iZ(a,1);
            Z(rz,cz) = X(dd,t);
            dd--;
            a--;
            
          }
        } 
        iF = inv_sympd(F.slice(t));
        Kt.slice(t) = T*Pt.slice(t)*Z.t()*iF;
        ut.col(t) = iF*v.col(t) - Kt.slice(t).t()*r.col(t);
        r.col(t-1) = Z.t()*ut.col(t) + T.t()*r.col(t);
        epshat.col(t) = H*ut.col(t);
        etahat.col(t) = Q1*R.t()*r.col(t);
        D.slice(t) = iF + Kt.slice(t).t()*N.slice(t)*Kt.slice(t);
        N.slice(t-1) =  Z.t()*D.slice(t)*Z + T.t()*N.slice(t)*T - Z.t()*Kt.slice(t).t()*N.slice(t)*T - T.t()*N.slice(t)*Kt.slice(t)*Z;
        V_eps.col(t) = (H - H*D.slice(t)*H).diag();
        V_eta.slice(t) = Q1 - Q1*R.t()*N.slice(t)*R*Q1; 
        
        dd = X.n_rows - 1;
        
        if (!(iR.n_rows == 2)) {
          while (iR(f,0) == 0){
            rz = iR(f,2);
            cz = iR(f,1);
            R(rz,cz) = X(dd,0);
            dd--;
            f--;
          }
        }
        
        if (!(iQ1.n_rows == 2)) {
          while (iQ1(e,0) == 0){
            rz = iQ1(e,2);
            cz = iQ1(e,1);
            Q1(rz,cz) = X(dd,0);
            dd--;
            e--;
          }
        }
        
        if (!(iH.n_rows == 2)) {
          while (iH(d,0) == 0){
            rz = iH(d,2);
            cz = iH(d,1);
            H(rz,cz) = X(dd,0);
            dd--;
            d--;
          }
        }
        
        if (!(iQ.n_rows == 2)) {
          while (iQ(c,0) == 0){
            rz = iQ(c,2);
            cz = iQ(c,1);
            Q(rz,cz) = X(dd,0);
            dd--;
            c--;
          }
        }
        
        if (!(iT.n_rows == 2)) {
          while (iT(b,0) == 0){
            rz = iT(b,2);
            cz = iT(b,1);
            T(rz,cz) = X(dd,0);
            dd--;
            b--;
          }
        }
        
        if (!(iZ.n_rows == 2)) {
          while (iZ(a,0) == 0){
            rz = iZ(a,2);
            cz = iZ(a,1);
            Z(rz,cz) = X(dd,0);
            dd--;
            a--;
          }
        } 
        
        iF = inv_sympd(F.slice(0));
        Kt.slice(0) = T*Pt.slice(0)*Z.t()*iF;
        ut.col(0) = iF*v.col(0) - Kt.slice(0).t()*r.col(0);
        epshat.col(0) = H*ut.col(0);
        etahat.col(0) = Q1*R.t()*r.col(0);
        D.slice(0) = iF + Kt.slice(0).t()*N.slice(0)*Kt.slice(0);
        V_eps.col(0) = (H - H*D.slice(0)*H).diag();
        V_eta.slice(0) = Q1 - Q1*R.t()*N.slice(0)*R*Q1;
      }
    }
  }
    if (Smo.size() == 2){
      if (((Rcpp::as<std::string>(Smo[0]) == "Smoother") and (Rcpp::as<std::string>(Smo[1]) == "Disturbance"))
            or ((Rcpp::as<std::string>(Smo[0]) == "Disturbance") and (Rcpp::as<std::string>(Smo[1]) == "Smoother"))){
        
        a = iZ.n_rows-2;
        b = iT.n_rows-2;
        c = iQ.n_rows-2;
        d = iH.n_rows-2;
        e = iQ1.n_rows-2;
        f = iR.n_rows-2;
        
        for (arma::uword t=n-1; t>=1; --t) {
          
          dd = X.n_rows - 1;
          
          if (!(iR.n_rows == 2)) {
            while (iR(f,0) == t){
              rz = iR(f,2);
              cz = iR(f,1);
              R(rz,cz) = X(dd,t);
              dd--;
              f--;
            }
          }
          
          if (!(iQ1.n_rows == 2)) {
            while (iQ1(e,0) == t){
              rz = iQ1(e,2);
              cz = iQ1(e,1);
              Q1(rz,cz) = X(dd,t);
              dd--;
              e--;
            }
          }
          
          if (!(iH.n_rows == 2)) {
            while (iH(d,0) == t){
              rz = iH(d,2);
              cz = iH(d,1);
              H(rz,cz) = X(dd,t);
              dd--;
              d--;
            }
          }
          
          if (!(iQ.n_rows == 2)) {
            while (iQ(c,0) == t){
              rz = iQ(c,2);
              cz = iQ(c,1);
              Q(rz,cz) = X(dd,t);
              dd--;
              c--;
            }
          }
          
          if (!(iT.n_rows == 2)) {
            while (iT(b,0) == t){
              rz = iT(b,2);
              cz = iT(b,1);
              T(rz,cz) = X(dd,t);
              dd--;
              b--;
            }
          }
          
          if (!(iZ.n_rows == 2)) {
            while (iZ(a,0) == t){
              rz = iZ(a,2);
              cz = iZ(a,1);
              Z(rz,cz) = X(dd,t);
              dd--;
              a--;
              
            }
          } 
          
          iF = inv_sympd(F.slice(t));
          Kt.slice(t) = T*Pt.slice(t)*Z.t()*iF;
          Lt.slice(t) = T - Kt.slice(t)*Z;
          ZiF = Z.t()*iF;
          r.col(t-1) = ZiF*v.col(t) + Lt.slice(t).t()*r.col(t); 
          N.slice(t-1) = ZiF*Z + Lt.slice(t).t()*N.slice(t)*Lt.slice(t);
          
          alphahat.col(t) = at.col(t) + Pt.slice(t)*r.col(t-1);   
          Vt.slice(t) = Pt.slice(t) - Pt.slice(t)*N.slice(t-1)*Pt.slice(t);
          
          ut.col(t) = iF*v.col(t) - Kt.slice(t).t()*r.col(t);
          epshat.col(t) = H*ut.col(t);
          etahat.col(t) = Q1*R.t()*r.col(t);
          D.slice(t) = iF + Kt.slice(t).t()*N.slice(t)*Kt.slice(t);
          V_eps.col(t) = (H - H*D.slice(t)*H).diag();
          V_eta.slice(t) = Q1 - Q1*R.t()*N.slice(t)*R*Q1;
        }
        
        dd = X.n_rows - 1;
        
        if (!(iR.n_rows == 2)) {
          while (iR(f,0) == 0){
            rz = iR(f,2);
            cz = iR(f,1);
            R(rz,cz) = X(dd,0);
            dd--;
            f--;
          }
        }
        
        if (!(iQ1.n_rows == 2)) {
          while (iQ1(e,0) == 0){
            rz = iQ1(e,2);
            cz = iQ1(e,1);
            Q1(rz,cz) = X(dd,0);
            dd--;
            e--;
          }
        }
        
        if (!(iH.n_rows == 2)) {
          while (iH(d,0) == 0){
            rz = iH(d,2);
            cz = iH(d,1);
            H(rz,cz) = X(dd,0);
            dd--;
            d--;
          }
        }
        
        if (!(iQ.n_rows == 2)) {
          while (iQ(c,0) == 0){
            rz = iQ(c,2);
            cz = iQ(c,1);
            Q(rz,cz) = X(dd,0);
            dd--;
            c--;
          }
        }
        
        if (!(iT.n_rows == 2)) {
          while (iT(b,0) == 0){
            rz = iT(b,2);
            cz = iT(b,1);
            T(rz,cz) = X(dd,0);
            dd--;
            b--;
          }
        }
        
        if (!(iZ.n_rows == 2)) {
          while (iZ(a,0) == 0){
            rz = iZ(a,2);
            cz = iZ(a,1);
            Z(rz,cz) = X(dd,0);
            dd--;
            a--;
          }
        } 
        
        iF = inv_sympd(F.slice(0));
        Kt.slice(0) = T*Pt.slice(0)*Z.t()*iF;
        Lt.slice(0) = T - Kt.slice(0)*Z;
        ZiF = Z.t()*iF;
        r0 = ZiF*v.col(0) + Lt.slice(0).t()*r.col(0); 
        N0 = ZiF*Z + Lt.slice(0).t()*N.slice(0)*Lt.slice(0);
        
        alphahat.col(0) = at.col(0) + Pt.slice(0)*r0;    
        Vt.slice(0) = Pt.slice(0) - Pt.slice(0)*N0*Pt.slice(0);
        
        ut.col(0) = iF*v.col(0) - Kt.slice(0).t()*r.col(0);
        epshat.col(0) = H*ut.col(0);
        etahat.col(0) = Q1*R.t()*r.col(0);
        D.slice(0) = iF + Kt.slice(0).t()*N.slice(0)*Kt.slice(0);
        
        V_eps.col(0) = (H - H*D.slice(0)*H).diag();
        V_eta.slice(0) = Q1 - Q1*R.t()*N.slice(0)*R*Q1;
      }
    }
    
  
  return List::create(Named("at") = at,
                      Named("Pt") = Pt,
                      Named("att") = att,
                      Named("Ptt") = Ptt,
                      Named("Yhat") = Ythat,
                      Named("K") = K,
                      Named("loglik") = -0.5*loglik,
                      Named("v") = v,
                      Named("F") = F,
                      Named("alphahat") = alphahat,
                      Named("Vt") = Vt,
                      Named("V_eps") = V_eps,
                      Named("V_eta") = V_eta,
                      Named("epshat") = epshat,
                      Named("etahat") = etahat  
                        
  );
}