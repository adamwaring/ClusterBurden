generate_input = function(dataset, type){

  x = dataset[order(dataset)]
  colnames(x) = c("p_pos", "ac", "an", "aff")
  x[,aff := ifelse(aff == unique(x$aff)[1], 1, 0)]
  ss = x[,floor(mean(an)/2),by=aff][[2]]
  x = x[rep(1:nrow(x), ac)]

  if(type == "pos_only"){

    input = x[,.(p_pos, aff)]

  } else if(type == "burden"){

    x = x[,.(p_pos, aff, car=1)]
    nvars = x[,sum(car), by=aff][[2]]

    n_0s = ss - nvars

    case_0s = matrix(rep(c(0, 1, 0), each=n_0s[1]), ncol=3)
    control_0s = matrix(rep(c(0, 0, 0), each=n_0s[2]), ncol=3)

    input = rbind(x[aff==1], case_0s, x[aff==0], control_0s, use.names=F)

  }

  return(input)

}


# generate gam model
make_model = function(input, type){

  if(type == "pos_only"){

    gam(aff ~ s(p_pos), method = "REML", family="binomial", data = input)

  } else if(type == "burden"){

    gam(aff ~ car + s(p_pos, by = car), method = "REML", family="binomial", data = input)

  }

}

# predict from model
predict_gam = function(mod, type, all_residues=F){

  input = data.table(mod$model)
  if(all_residues){
    input = data.table(p_pos=1:max(input$p_pos), aff=1, car=1)
  } else if(type == "pos_only"){
    input = input[aff==1]
  } else if(type == "burden"){
    input = input[car==1 && aff == 1]
  }

  or = rowSums(predict(mod, newdata=input, type="terms"))
  se = predict(mod, newdata=input, se.fit=T)$se.fit
  upr <- or + (2 * se)
  lwr <- or - (2 * se)
  or = exp(or); upr = exp(upr); lwr = exp(lwr)

  input[, c("or", "upr", "lwr") := .(or, upr, lwr)]

}

# plot model predictions
plot_gam = function(preds, all_residues){

  g1 = ggplot(preds, aes(x=p_pos, y = or)) + geom_point(size=2) +
    labs(x = "Residue position in linear protein sequence", y = "Predicted odds-ratio") +
    theme_bw(20) + theme(panel.border = element_blank(),
                         axis.line = element_line(),
                         plot.margin = unit(rep(0.001, 4), "cm"),
                         panel.grid = element_blank())

  if(all_residues){
    g1 = g1 + geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.1)
  } else {
    g1 = g1 + geom_errorbar(aes(ymin=lwr, ymax=upr))
  }

  return(g1)

}

gam_pipeline = function(dataset, type="pos_only", all_residues=F){

  input = generate_input(dataset, type)
  mod = make_model(input, type)
  preds = predict_gam(mod, type, all_residues)
  plot_gam(preds, all_residues)

}
