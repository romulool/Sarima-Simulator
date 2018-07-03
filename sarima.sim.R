setwd("/media/romulo/Windows/Agenda/UFPR/Estatística/Laboratórios/codigo")

sarima.sim = function(n = 1,
                      ar = NULL,
                      sar = NULL,
                      ma = NULL,
                      sma = NULL,
                      D = 0,
                      d = 0,
                      s = 1,
                      sigma2 = 1,
                      delta = 0){

# defaults
#n = 1
#ar = NULL
#sar = NULL
#ma = NULL
#sma = NULL
#D = 0
#d = 0
#s = 1
#sigma2 = 1
#delta = 0

# entradas
#n = 60
#ar = c(-0.2962)
#ar = c(-0.9)
#sar = c(0.82,0.26)
#sar = c(0.82)
#ma = c(0.76,0.49)
#ma = c(0.2)
#sma = c(-1.9306,1.3063,-0.2875)
#sma = c(-0.2875)
#D = 3
#d = 1
#s = 4
#sigma2 = 160
#delta = 300

# Transformação para formato de lista dos parâmetros ar, sar, ma e sma
l = list(c(1,0))
if (!(is.null(ar))){
  for (i in 1:length(ar)){
    l = append(l,list(c(-ar[i],i)))
  }
}
ar = l

l = list(c(1,0))
if (!(is.null(sar))){
  for (i in 1:length(sar)){
    l = append(l,list(c(-sar[i],i*s)))
  }
}
sar = l

l = list(c(1,0))
if (!(is.null(ma))){
  for (i in 1:length(ma)){
    l = append(l,list(c(-ma[i],i)))
  }
}
ma = l

l = list(c(1,0))
if (!(is.null(sma))){
  for (i in 1:length(sma)){
    l = append(l,list(c(-sma[i],i*s)))
  }
}
sma = l

l = list(c(1,0))
if (d > 0){
  for (i in 1:d){
    l = append(l,list(c(((-1)^i)*((factorial(d))/(factorial(i)*factorial(d - i))),i)))
  }
}
reg_dif = l

l = list(c(1,0))
if (D > 0){
  for (i in 1:D){
    l = append(l,list(c(((-1)^i)*((factorial(D))/(factorial(i)*factorial(D - i))),i*s)))
  }
}
saz_dif = l

# 1st operação: produto de polinômios

ar
sar
ar_sar = list()
for(i in 1:length(ar)){
  for(j in 1:length(sar)){
    ar_sar = append(ar_sar,list(c(ar[[i]][1] * sar[[j]][1], ar[[i]][2] + sar[[j]][2])))
  }
}
ar_sar

ar_sar
saz_dif
ar_sar_D = list()
for(i in 1:length(ar_sar)){
  for(j in 1:length(saz_dif)){
    ar_sar_D = append(ar_sar_D,list(c(ar_sar[[i]][1] * saz_dif[[j]][1], ar_sar[[i]][2] + saz_dif[[j]][2])))
  }
}
ar_sar_D

ar_sar_D
reg_dif
ar_sar_D_d = list()
for(i in 1:length(ar_sar_D)){
  for(j in 1:length(reg_dif)){
    ar_sar_D_d = append(ar_sar_D_d,list(c(ar_sar_D[[i]][1] * reg_dif[[j]][1], ar_sar_D[[i]][2] + reg_dif[[j]][2])))
  }
}
ar_sar_D_d

ma
sma
ma_sma = list()
for(i in 1:length(ma)){
  for(j in 1:length(sma)){
    ma_sma = append(ma_sma,list(c(ma[[i]][1] * sma[[j]][1], ma[[i]][2] + sma[[j]][2])))
  }
}
ma_sma

# 2nd operação: soma associativa
v = NULL
for (i in 1:length(ar_sar_D_d)){
  v = append(v,ar_sar_D_d[[i]][2])  
}
values_z = as.integer(names(table(v)))

lozl = list()
k = 1
for (i in values_z){
  soma = 0
  for (j in 1:length(ar_sar_D_d)){
    if (ar_sar_D_d[[j]][2] == as.numeric(i)){
      soma = soma + ar_sar_D_d[[j]][1]
    }
  }
  lozl = append(lozl,list(c(soma,i)))
  k = k + 1
}
lozl

v = NULL
for (i in 1:length(ma_sma)){
  v = append(v,ma_sma[[i]][2])  
}
values_e = as.integer(names(table(v)))

loe = list()
k = 1
for (i in values_e){
  soma = 0
  for (j in 1:length(ma_sma)){
    if (ma_sma[[j]][2] == as.numeric(i)){
      soma = soma + ma_sma[[j]][1]
    }
  }
  loe = append(loe,list(c(soma,i)))
  k = k + 1
}
loe

# 3rd operação: isolar z_t à esquerda da equação
lozr = list()
if (length(lozl) > 1){
  for (i in 2:length(lozl)){
    lozr = append(lozr,list(lozl[[i]]))
    lozr[[i - 1]][1] = -lozr[[i - 1]][1]
  }
}
lozr

# 4th operação: gerar valores aleatórios
rand = rnorm(n, mean = 0, sd = sqrt(sigma2))
limit = max(max(values_z),max(values_e)) + 1
serie = rep(delta,limit)
error = rep(0,limit)

for (i in 1:n){
  z = delta
  if (length(lozr) > 0){  # verificando se há parte ar
    for (j in 1:length(lozr)){
      z = z + lozr[[j]][1]*serie[length(serie) - as.integer(lozr[[j]][2]) + 1]
    }
  }
  for (j in 1:length(loe)){
    if (as.integer(loe[[j]][2]) == 0){
      z = z + loe[[j]][1]*rand[i]
    }else{
      z = z + loe[[j]][1]*error[length(error) - as.integer(loe[[j]][2])]
    }
  }
  serie = append(serie,z)
  error = append(error,rand[i])
}

serie = serie[seq(limit + 1,length(serie))]
error = error[seq(limit + 1,length(error))]

out = list(ts(serie),ts(error))
names(out) = c("serie","error")

return(out)
}


sim = sarima.sim(n = 120,
                 ar = -0.2,
                 ma = 0.8,
                 sma = c(-1.93),
                 sar = c(0.2),
                 D = ,
                 s = 4,
                 sigma2 = 25
                 #delta = 0
                 )
plot(sim$serie)
#lines(1:n,sim$serie)
plot(sim$error)
#lines(1:n,sim$error)

# Ajuste de verificação
fit = arima(sim$serie,order = c(0,0,0),seasonal = list(order=c(0,0,0)));fit
#fit = arima(serie,order = c(1,0,0));fit
acf(sim$serie)
pacf(sim$serie)
