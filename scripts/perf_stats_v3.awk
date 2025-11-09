BEGIN { FS=","; saw_time=0 }
function isnum(x){ return (x ~ /^-?[0-9]+(\.[0-9]+)?$/) }
function sort2(A,B,n, a,b,v0,v1){ for(a=2;a<=n;a++){ v0=A[a]; v1=B[a]; b=a-1; while(b>=1&&A[b]>v0){A[b+1]=A[b]; B[b+1]=B[b]; b--} A[b+1]=v0; B[b+1]=v1 } }
function pidx(n,p, x){ x=int((p/100.0)*n+0.999999); if(x<1)x=1; if(x>n)x=n; return x }
function stats(lbl,X,n, k,sum,s2v,mean,sd,minv,med,p90,p99){
  if(n==0){ printf("%s: none\n",lbl); return }
  for(k=1;k<=n;k++) T[k]=X[k]
  sum=s2v=0; minv=T[1]
  for(k=1;k<=n;k++){ sum+=T[k]; s2v+=T[k]*T[k]; if(T[k]<minv) minv=T[k] }
  sort2(T,T,n); mean=sum/n; sd=(n>1)?sqrt(s2v/n-mean*mean):0
  med=T[pidx(n,50)]; p90=T[pidx(n,90)]; p99=T[pidx(n,99)]
  printf("%s: mean=%.6f sd=%.6f min=%.6f median=%.6f p90=%.6f p99=%.6f\n", lbl, mean, sd, minv, med, p90, p99)
}
$0 ~ /,/ {
  ev = ($4 != "" ? $4 : $3); cnt = $2; ts=$1
  if (ts ~ /^[0-9.]+$/) {
    if (ts != cur_ts) {
      cur_ts = ts
      cur_idx = i + 1
    }
    saw_time=1
  }
  if ((ev ~ /mem_inst_retired\.all_stores/) && isnum(cnt)) {
    if (cur_idx == 0) cur_idx = i + 1
    i = cur_idx
    S[i]=cnt+0
    next
  }
  if (!isnum(cnt)) next
  target = (cur_idx > 0 ? cur_idx : i)
  if (target <= 0) next
  if (((ev ~ /cpu_core\/cycles\//) || (ev ~ /^cycles$/)))                          CY[target]=cnt+0
  else if (((ev ~ /cpu_core\/instructions\//) || (ev ~ /^instructions$/)))         IN[target]=cnt+0
  else if (ev ~ /resource_stalls\.sb/)                                             R[target]=cnt+0
  else if (ev ~ /ld_blocks\.store_forward/)                                        F[target]=cnt+0
  else if (ev ~ /cache-references/)                                                C[target]=cnt+0
  else if (ev ~ /cache-misses/)                                                    M[target]=cnt+0
}
END{
  if(i==0){ print "no rows found"; exit }
  warmup = (ENVIRON["WARMUP"] ~ /^[0-9]+$/) ? ENVIRON["WARMUP"]+0 : (saw_time ? 8 : 0)
  trim   = (ENVIRON["TRIM"]   ~ /^[0-9]+$/) ? ENVIRON["TRIM"]+0   : (saw_time ? 2 : 0)
  n=0
  for(k=warmup+1;k<=i;k++){
    if (S[k]>0 && C[k]>0){
      sps[++n] = R[k]/S[k]
      mr[n]    = M[k]/C[k]
      if (k in CY) cps[n] = CY[k]/S[k]
      if (k in IN) ips[n] = IN[k]/S[k]
      if (k in F)  sfb[n] = F[k]/S[k]
      volS[n]=S[k]; volC[n]=C[k]
    }
  }
  if(n==0){ print "no usable rows after warmup"; exit }
  for(k=1;k<=n;k++){ idx[k]=k; key[k]=sps[k] } sort2(key,idx,n)
  from=1+trim; to=n-trim; if(from>to){from=1;to=n}
  used=0
  for(t=from;t<=to;t++){
    k=idx[t]; used++
    fsps[used]=sps[k]; fmr[used]=mr[k]; fS[used]=volS[k]; fC[used]=volC[k]
    if (k in cps) fcps[used]=cps[k]
    if (k in ips) fips[used]=ips[k]
    if (k in sfb) fsfb[used]=sfb[k]
  }
  printf("rows_total=%d warmup_skipped=%d rows_used=%d trimmed_each_side=%d autodetect_interval=%s\n",
         i, warmup, used, trim, saw_time ? "yes" : "no")
  sumS=sumC=0; for(k=1;k<=used;k++){ sumS+=fS[k]; sumC+=fC[k] }
  printf("avg_stores_per_row=%.3f  avg_refs_per_row=%.3f\n", sumS/used, sumC/used)
  stats("stalls_per_store", fsps, used)
  stats("miss_rate",        fmr,  used)
  have_cps=have_ips=have_sfb=0
  for(k=1;k<=used;k++){ if(k in fcps) have_cps=1; if(k in fips) have_ips=1; if(k in fsfb) have_sfb=1 }
  if(have_cps){ j=0; for(k=1;k<=used;k++) if(k in fcps) TMP[++j]=fcps[k]; stats("cycles_per_store", TMP, j) } else print "cycles_per_store: cycles not present in file"
  if(have_ips){ j=0; for(k=1;k<=used;k++) if(k in fips) TMP[++j]=fips[k]; stats("instructions_per_store", TMP, j) } else print "instructions_per_store: instructions not present in file"
  if(have_sfb){ j=0; for(k=1;k<=used;k++) if(k in fsfb) TMP[++j]=fsfb[k]; stats("store_forward_blocks_per_store", TMP, j) } else print "store_forward_blocks_per_store: not present in file"
  tsumS=tsumC=tsumSt=tsumMr=tsumCPS=tsumIPS=tsumSFB=0
  for(k=1;k<=used;k++){
    tsumS+=fS[k]; tsumC+=fC[k]
    tsumSt+=fsps[k]*fS[k]
    tsumMr+=fmr[k]*fC[k]
    if(k in fcps) tsumCPS+=fcps[k]*fS[k]
    if(k in fips) tsumIPS+=fips[k]*fS[k]
    if(k in fsfb) tsumSFB+=fsfb[k]*fS[k]
  }
  printf("ratio_of_sums: stalls_per_store=%.6f  miss_rate=%.6f", (tsumS>0?tsumSt/tsumS:0), (tsumC>0?tsumMr/tsumC:0))
  if(have_cps) printf("  cycles_per_store=%.6f", (tsumS>0?tsumCPS/tsumS:0))
  if(have_ips) printf("  instructions_per_store=%.6f", (tsumS>0?tsumIPS/tsumS:0))
  if(have_sfb) printf("  store_forward_blocks_per_store=%.6f", (tsumS>0?tsumSFB/tsumS:0))
  printf("\n")
}
