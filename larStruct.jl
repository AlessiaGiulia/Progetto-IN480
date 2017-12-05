using PyCall
@pyimport scipy


function prepKey (args)
	return  
end

function fixedPrec(PRECISION)
	function fixedPrec0(value)
		out=round(value,PRECISION)
		if out==-0.0
			out=0.0
		end
		return string(out)
	end
	return fixedPrec0
end

 

function vcode(PRECISION=4)

            function vcode0(vect)

                        Return prepKey(map(fixedPrec(PRECISION),vect))

            End

            Return vcode0

end

function t(args…)

            d=length(args)

            mat=eye(d+1)

            for k in range(1,d)

                        mat[k,d+1]=args[k]

            end

            return mat

end

 

function s(args…)

            d=length(args)

            mat=eye(d+1)

            for k in range(1,d)

                        mat[k,k]=args[k]

            end

            return mat

end

function removeDups(CW)
	CW=[(map(tuple,CW))]
	return CW
end 

 
function larRemoveVertices(V,FV)
	vertDict= Dict()
	index,defaultValue,CW,W,FW = 0,-1,[],[],[]
	for (k,incell) in enumerate(FV)
		outcell=[]
		for v in incell
			key=vcode(4)(V[v])
			if get(vertDict,key,defaultValue)==defaultValue
				index =index+1
				vertDict[key]=index
				outcell =append!(outcell,index])
				W= append!(W,eval(key))                   
			else
				outcell= append!(outcell,vertDict[key])
			end
		end
		FW =append!(FW,[outcell])
	end
	return W,FW
end

function larBoundary(self)
	data=struct2lar(self)
	if length(data)==3
		V,FV,EV=data
		return V,FV,EV
	else
	return "<Struct name $self.__name__() : boundary non computable"   #controllare che effettivamente ritorna il nome della struct
	end
end



 

function struct2lar(structure,metric=ID)

            listOfModels=evalStruct(structure)

            vertDict= Dict()

            index,defaultValue,CW,W,FW = 0,-1,[],[],[]

            for model in listOfModels

                        if  length(model)==2

                                   V,FV=model

                        elseif

                                   V,FV,EV=model

                        end

                        for (k,incell) in enumerate(FV)

                                   outcell=[]

                                   for v in incell

                                   	key=vcode(4)(V[v])

                                   	if get(vertDict,key,defaultValue)==defaultValue

                                        	indexx =index+1

                                                vertDict[key]=index

                                                outcell =append!(outcell,index])
						W= append!(W,eval(key))                   
                                           else

                                                 outcell= append!(outcell,vertDict[key])

                                               end

                                   CW =append!(CW,[outcell])

                                   End

                        end

                        If length(model)==3

                                   outcell=[]

                                   for v in incell

                                               key=vcode(4)(V[v])

                                               if get(vertDict,key,defaultValue)==defaultValue

                                                           index =+#vedere se devo riscrivere index

                                                           vertDict[key]=index

                                                           outcell =append!(outcell,[index])#se da errore vedere un

W= append!(W,eval(key))                   metodo per le liste

                                               else

                                                           outcell= append!(outcell,vertDict[key])

                                               end

                                               FW =append!(FW,[outcell])

                                   End

                        end

                        if length(model)==2

                                   elseif length(CW[0])==2

                                               CW= list(set(AA(tuple)(AA(sorted)(CW))))

                                    Else

                                               CW= removeDups(CW)

                                   Return metric(W),CW

End

If length(model)==3

            FW= list(set(AA(tuple)(AA(sorted)(FW))))

                                    CW= removeDups(CW)

                                   Return metric(W),CW,FW

                        end

            end

end

 

 

 


=======
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

>>>>>>> origin/master

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

