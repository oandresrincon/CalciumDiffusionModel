
# coding: utf-8

# In[3]:


def do_timestep(u0, u, D,dl,dr,dt,m,n):

    # first row
    
    u_down = 0#D * ( u0[1,0] - u0[0,0] ) / dr#0
    u_up =  0
    u_right =  D * ( u0[0,1] - u0[0,0] ) / dl 
    u_left = 0
    u[0,0] = u0[0,0] + ( dt * ( u_down + u_up + u_right + u_left ) ) 

    u_down = 0#D * ( u0[1,1:m-1] - u0[0,1:m-1] ) / dr#0  
    u_up =  0
    u_right =  D * ( u0[0,2:m] - u0[0,1:m-1] ) / dl 
    u_left = D * ( u0[0,0:m-2] - u0[0,1:m-1] ) / dl
    u[0,1:m-1] = u0[0,1:m-1] + ( dt * ( u_down + u_up + u_right + u_left ) ) 
            
    u_down =  0#D * ( u0[1,m-1] - u0[0,m-1] ) / dr#0
    u_up =  0
    u_right = 0 
    u_left = D * ( u0[0,m-2] - u0[0,m-1] ) / dl
    u[0,m-1] = u0[0,m-1] + ( dt * ( u_down + u_up + u_right + u_left ) )    
    
    # second row
    
    u_down = D * ( u0[2,0] - u0[1,0] ) / dr  
    u_up = 0  
    u_right =  D * ( u0[1,1] - u0[1,0] ) / dl 
    u_left = 0 
    u[1,0] = u0[1,0] + ( dt * ( u_down + u_up + u_right + u_left ) ) 

            
    u_down = D * ( u0[2,1:m-1] - u[1,1:m-1] ) / dr 
    u_up = 0  
    u_right = D * ( u0[1,2:m] - u[1,1:m-1] ) / dl 
    u_left = D * ( u0[1,0:m-2] - u[1,1:m-1] ) / dl
    u[1,1:m-1] = u0[1,1:m-1] + ( dt * ( u_down + u_up + u_right + u_left ) ) 
                    
            
    u_down = D * ( u0[2,m-1] - u0[1,m-1] ) / dr 
    u_up = 0#D * ( u0[i-1,j] - u0[i,j] ) / dr 
    u_right = 0#D * ( u0[i,j+1] - u0[i,j] ) / dl 
    u_left = D * ( u0[1,m-2] - u0[1,m-1] ) / dl
    u[1,m-1] = u0[1,m-1] + ( dt * ( u_down + u_up + u_right + u_left ) )     
    
    # rigth column
    
    u_down = D * ( u0[3:n,m-1] - u0[2:n-1,m-1] ) / dr  
    u_up = D * ( u0[1:n-2,m-1] - u0[2:n-1,m-1] ) / dr 
    u_right = 0 
    u_left = D * ( u0[2:n-1,m-2] - u0[2:n-1,m-1] ) / dl
    u[2:n-1,m-1] = u0[2:n-1,m-1] + ( dt * ( u_down + u_up + u_right + u_left ) )
    
    # left row
    
    u_down = D * ( u0[3:n,0] - u0[2:n-1,0] ) / dr 
    u_up = D * ( u0[1:n-2,0] - u0[2:n-1,0] ) / dr 
    u_right = D * ( u0[2:n-1,1] - u0[2:n-1,0] ) / dl 
    u_left = 0
    u[2:n-1,0] = u0[2:n-1,0] + ( dt * ( u_down + u_up + u_right + u_left ) )         
    
    #down row
    
    u_down = 0
    u_up = D * ( u0[n-2,m-1] - u0[n-1,m-1] ) / dr 
    u_right = 0 
    u_left = D * ( u0[n-1,m-2] - u0[n-1,m-1] ) / dl
    u[n-1,m-1] = u0[n-1,m-1] + ( dt * ( u_down + u_up + u_right + u_left ) ) 
            
    u_down = 0
    u_up = D * ( u0[n-2,1:m-1] - u0[n-1,1:m-1] ) / dr 
    u_right =  D * ( u0[n-1,2:m] - u0[n-1,1:m-1] ) / dl
    u_left = D * ( u0[n-1,0:m-2] - u0[n-1,1:m-1] ) / dl
    u[n-1,1:m-1] = u0[n-1,1:m-1] + ( dt * ( u_down + u_up + u_right + u_left ) ) 
                                
    u_down = 0
    u_up = D * ( u0[n-2,0] - u0[n-1,0] ) / dr 
    u_right =  D * ( u0[n-1,1] - u0[n-1,0] ) / dl
    u_left = 0
    u[n-1,0] = u0[n-1,0] + ( dt * ( u_down + u_up + u_right + u_left ) )            
    
    #inside
    
    u_down = D * ( u0[3:n,1:m-1] - u0[2:n-1,1:m-1] ) / dr  
    u_up = D * ( u0[1:n-2,1:m-1] - u0[2:n-1,1:m-1] ) / dr 
    u_right = D * ( u0[2:n-1,2:m] - u0[2:n-1,1:m-1] ) / dl 
    u_left = D * ( u0[2:n-1,0:m-2] - u0[2:n-1,1:m-1] ) / dl
    u[2:n-1,1:m-1] = u0[2:n-1,1:m-1] + ( dt * ( u_down + u_up + u_right + u_left ) ) 

    u0 = u.copy()
    return u0, u

