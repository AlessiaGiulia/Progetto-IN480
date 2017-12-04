using PyCall
@pyimport scipy



--------------------------------------------------------------------------------------------------------------------------------------------

larApply (affineMatrix)
  larApply0(model)
    if length(model)==2
      V,CV=model
    elseif length(model)==3
      V,CV,FV = model
    end
    V=reshape(for v in V append!(v,1.0))*(transpose(affineMatrix))
  return larApply0


  

