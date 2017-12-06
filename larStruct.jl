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
		return prepKey(map(fixedPrec(PRECISION),vect))
	end
	return vcode0
end


function t(args…)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
        	mat[k,d+1]=args[k]
	end
	return mat.view(Mat) #vedere view dopo aver fatto la classe Mat
end

 

function s(args…)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
		mat[k,k]=args[k]
	end
	return mat
end


function larEmbed(k)
	function larEmbed0(model)
		if k>0
			model[1]=[append!(v,fill(0.0,k)) for v in model[1]]
		elseif k<0
			for n in range(1,length(model[1][1])+k)	#vedere se esiste una funzione più pratica
				for v in model[1]
					pop!(v)
				end
			end
		end3
		return model
	end
	return larEmbed0
end


function removeDups(CW)
	CW=collect(tuple(CW...))
	CWs=collect(tuple(sort(CW,rev=true)...))
	no_duplicates=Dict()
	for f in CWs
		no_duplicates[f] = []
	end
	for f in CW
	no_duplicates[tuple(sort(collect(f))...)]=[f]
	end
	CW=[f[1] for f in values(no_duplicates)]
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
		elseif lenght(model)==3
			V,FV,EV=model
		end
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
			CW =append!(CW,[outcell])
		end
		if length(model)==3
			for (k,incell) in enumerate(FV)
				outcell=[]
				for v in incell
					key=vcode(4)(V[v])
					if get(vertDict,key,defaultValue)==defaultValue
						index =index+1
						vertDict[key]=index
						outcell =append!(outcell,[index])
						W=append!(W,eval(key))                   
					else
						outcell= append!(outcell,vertDict[key])
					end
				end
				FW =append!(FW,[outcell])
			end
		end
	end
	if length(model)==2
		elseif length(CW[0])==2
			CW=list(set(AA(tuple)(AA(sorted)(CW))))#da cambiare
		else
			CW=removeDups(CW)
		return metric(W),CW
	end
	if length(model)==3
		FW=list(set(AA(tuple)(AA(sorted)(FW))))#da cambiare
		CW=removeDups(CW)
		return metric(W),CW,FW
	end
end

function embedTraversal(cloned,obj,n,suffix)
	for i in range(1,length(obj))
		if (isa(obj[i],tuple) ||isa(obj[i],Array))&& length(obj[i])==2 #potrebbe dare problemi in ndarray
			V,EV=obj[i]
			#V=[append!(v,fill(0.0,n)) for v in V]    #provare a vedere se girano entrambe allo stesso modo
			dimadd=fill([0.0],n)	#controllare se esiste una funzione
			for k in dimadd
				for v in V
					append!(v,k)
				end
			end
			append!(cloned.body,(V,EV))	#controllare dopo aver fatto la classe Struct
		elseif (isa(obj[i],tuple) ||isa(obj[i],Array))&& length(obj[i])==3 #potrebbe dare problemi in ndarray
			V,FV,EV=obj[i]
			dimadd=fill([0.0],n)	#controllare se esiste una funzione
			for k in dimadd
				for v in V
					append!(v,k)
				end
			end
			append!(cloned.body,(V,FV,EV))	#controllare dopo aver fatto la classe Struct
		elseif isa(obj[i],Mat)
			mat=obj[i]
			d,d=size(mat)
			newMat=eye(d+n*1)
			for h in range(1,d-1)
				for k in range(1,d-1)
					newMat[h,k]=mat[h,k]
				end
				newMat[h,d-1+n*1]=mat[h,d-1]
			end
			append!(cloned.body,newMat.view(Mat))#controllare view 
		elseif isa(obj[i],Struct)	#controllare dopo aver fatto la classe Struct
			newObj=Struct()
			newObj.box=hstack((obj[i].box,[fill([0],n),fill([0],n)]))	#vedere a chi corrisponde hstack
			newObj.category=obj[i].category
			append!(cloned.body,embedTraversal(newObj,obj[i],n,suffix))	#controllare dopo aver fatto la classe Struct
		end
	end
	return cloned
end

function embedStruct(n)
	function embedStruct0(struct,suffix="New")
		if n==0
			return struct, length(struct.box[1])
		end
		cloned=Struct()
		cloned.box=hstack((struct.box,[fill([0],n),fill([0],n)]))	#vedere a chi corrisponde hstack
		cloned.name=string(object_id(cloned))
		cloned.category=struct.category
		cloned.dim=struct.dim+n
		cloned=embedTraversal(cloned,struct,n,suffix)
		return cloned
	end
end

 

 

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
    return [v[:-1] for v in V],CV		#una cosa simile la avevo anche io...ho fatto un for(in larEmbed),ricordiamoci di confrontare!

  elseif len(model)==3
    return [v[:-1] for v in V],CV,FV
 end
end 


	
