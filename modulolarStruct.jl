using PyCall
@pyimport scipy



--------------------------------------------------------------------------------------------------------------------------------------------

function larApply (affineMatrix)
  larApply0(model)
    if length(model)==2
      V,CV=model
    elseif length(model)==3
      V,CV,FV = model
    end
    for v in V 
      append!(v,1.0) 
    end   
    V=(vec(v))*(transpose(affineMatrix))
  return larApply0

  if len(model)==2
    return [v[:-1] for v in V],CV

  elseif len(model)==3
    return [v[:-1] for v in V],CV,FV
 end
end 

