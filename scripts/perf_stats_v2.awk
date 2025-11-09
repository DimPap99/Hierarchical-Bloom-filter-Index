BEGIN {
  FS=","; saw_time=0
}
function isnum(x){ return (x ~ /^-?[0-9]+(\.[0-9]+)?$/) }
function is_allstores(ev){ return (ev ~ /(cpu_(core|atom)\/)?mem_inst_retired\.all_stores\/?/) }
function is_stalls_sb(ev){ return (ev ~ /(cpu_(core|atom)\/)?resource_stalls\.sb\/?/) }
function is_cache_refs(ev){ return (ev ~ /(cpu_(core|atom)\/)?cache-references\/?/) }
function is_cache_miss(ev){ return (ev ~ /(cpu_(core|atom)\/)?cache-misses\/?/) }
function is_store_fwd(ev){ return (ev ~ /(cpu_(core|atom)\/)?ld_blocks\.store_forward\/?/) }
function sort2(A, B, n,    a,b,v0,v1){
  for (a=2;a<=n;a++){ v0=A[a]; v1=B[a]; b=a-1; while(b>=1 && A[b]>v0){A[b+1]=A[b]; B[b+1]=B[b]; b--} A[b+1]=v0; B[b+1]=v1 }
}
function pidx(n,p,    x){ x=int((p/100.0)*n+0.999999); if(x<1)x=1; if(x>n)x=n; return x }
function stats(label, X, n,    k,sum,sumsq,mean,sd,minv,med,p90,p99){
  if(n==0){ printf("%s: none\n",label); return }
  for(k=1;k<=n;k++) T[k]=X[k]
  sum=sumsq=0; minv=T[1]
  for(k=1;k<=n;k++){ sum+=T[k]; sumsq+=T[k]*T[k]; if(T[k]<minv) minv=T[k] }
  sort2(T,T,n)
  mean=sum/n; sd=(n>1)? sqrt(sumsq/n-mean*mean):0
  med=T[pidx(n,50)]; p90=T[pidx(n,90)]; p99=T[pidx(n,99)]
  printf("%s: mean=%.6f sd=%.6f min=%.6f median=%.6f p90=%.6f p99=%.6f\n", label, mean, sd, minv, med, p90, p99)
}
$0 ~ /,/ {
  ev=$4; if (ev=="") ev=$3
  cnt=$2
  if (is_allstores(ev) && isnum(cnt)) { ++i; S[i]=cnt+0; if ($1 ~ /^[0-9.]+$/) saw_time=1; next }
  if (i==0) next
  if (is_stalls_sb(ev) && isnum(cnt)) R[i]=cnt+0
  else if (is_cache_refs(ev) && isnum(cnt)) C[i]=cnt+0
  else if (is_cache_miss(ev) && isnum(cnt)) M[i]=cnt+0
  else if (is_store_fwd(ev) && isnum(cnt)) F[i]=cnt+0
}
END{
  if(i==0){ print "no rows found"; exit }
  warmup = (ENVIRON["WARMUP"] ~ /^[0-9]+$/) ? ENVIRON["WARMUP"]+0 : (saw_time ? 8 : 0)
  trim   = (ENVIRON["TRIM"]   ~ /^[0-9]+$/) ? ENVIRON["TRIM"]+0   : (saw_time ? 2 : 0)

  n=0
  for(k=warmup+1;k<=i;k++){
    if (S[k]>0 && C[k]>0){
      sps[++n]=R[k]/S[k]
      mr[n]=M[k]/C[k]
      if(F[k]!="") sfb[n]=F[k]/S[k]
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
    if (k in sfb) { fsfb[used]=sfb[k]; haveSFB=1 }
  }

  printf("rows_total=%d warmup_skipped=%d rows_used=%d trimmed_each_side=%d autodetect_interval=%s\n",
         i, warmup, used, trim, saw_time ? "yes" : "no")

  sumS=sumC=0; for(k=1;k<=used;k++){ sumS+=fS[k]; sumC+=fC[k] }
  printf("avg_stores_per_row=%.3f  avg_refs_per_row=%.3f\n", sumS/used, sumC/used)

  stats("stalls_per_store", fsps, used)
  stats("miss_rate",        fmr,  used)

  if (haveSFB) { j=0; for(k=1;k<=used;k++) if (k in fsfb || fsfb[k]!="") TMP[++j]=fsfb[k]; stats("store_forward_blocks_per_store", TMP, j) }
  else { print "store_forward_blocks_per_store: not present in file" }

  tsumS=tsumC=tsumSt=tsumMr=0
  for(k=1;k<=used;k++){ tsumS+=fS[k]; tsumC+=fC[k]; tsumSt+=fsps[k]*fS[k]; tsumMr+=fmr[k]*fC[k] }
  printf("ratio_of_sums: stalls_per_store=%.6f  miss_rate=%.6f\n",
         (tsumS>0?tsumSt/tsumS:0), (tsumC>0?tsumMr/tsumC:0))
}
