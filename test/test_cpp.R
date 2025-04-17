library(lme4)
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
fit <- lmer(mform, data = Pixel)
precomp <- limestest:::get_precomp_lmer(fit)
H <- do.call(cbind, precomp$Hlist)
psi_hat <- limestest:::get_psi_hat_lmer(fit)

Psi_cpp <- limestest:::Psi_from_H_cpp(psi_mr = psi_hat[-length(psi_hat)], H = H)
Psi_R <- limestest:::Psi_from_Hlist(psi_mr = psi_hat[-length(psi_hat)], Hlist = precomp$Hlist)

cat(max(abs(Psi_cpp - Psi_R)), "\n")
