
function my_lens, theta, pars, z_weight=z_weight, alpha_only=alpha_only
; This is where we put the specific lens model that will be
; used. We can also combine different analytical models to create the
; full lens model.

  lens=theta*0d

; For now we just have a set of NIEs that all contribute to the lens model

  nnie=n_elements(pars)/8L ; How many nie models?
  if nnie gt 0 then lens = nie_lens(theta,pars[0:7],$
                                    z_weight=z_weight,alpha_only=keyword_set(alpha_only))

  for i=1,nnie-1 do lens += nie_lens(theta,pars[i*8L:(i+1)*8L - 1L],$
                                     z_weight=z_weight,alpha_only=keyword_set(alpha_only))

  return, lens
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;








; older version...

;function my_lens, theta, pars, z_weight=z_weight, alpha_only=alpha_only
; This is where we put the specific lens model that will be
; used. We can also combine different analytical models to create the
; full lens model.

;  if keyword_set(alpha_only) then a=1 else a=0

;  nie0_pars=pars[0:7]
;  nie1_pars=pars[8:23]
;  nie2_pars=pars[24:31]


;  lens=nie_lens(theta,nie0_pars,z_weight=z_weight,alpha_only=a)+$
;       nie_lens(theta,nie1_pars,z_weight=z_weight,alpha_only=a);+$
;       nie_lens(theta,nie2_pars,z_weight=z_weight,alpha_only=a)

;  return, lens
;end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


