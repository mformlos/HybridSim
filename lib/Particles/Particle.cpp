#include "Particle.h" 

Particle::Particle() : 
    Position { Vector3d::Zero() }, 
    Velocity { Vector3d::Zero() }, 
    Mass { 1.0 } { }
    
Particle::Particle(double aMass) : 
    Position { Vector3d::Zero() }, 
    Velocity { Vector3d::Zero() }, 
    Mass { aMass } { }
    
Particle::Particle(Vector3d aPosition, Vector3d aVelocity, double aMass) :
	Position { aPosition },
	Velocity { aVelocity },
	Mass { aMass } { }

MPCParticle::MPCParticle() : 
    Particle() { }
    
MPCParticle::MPCParticle(double aMass) : 
    Particle(aMass) { }
    
MPCParticle::MPCParticle(Vector3d aPosition, Vector3d aVelocity, double aMass) :
    Particle(aPosition, aVelocity, aMass) { }
    

	

MDParticle::MDParticle() : 
    MPCParticle(), 
    Force { Vector3d::Zero() }, 
    VerletPosition { Vector3d::Zero() }, 
    VerletList { } { }


MDParticle::MDParticle(double aMass) : 
    MPCParticle(aMass), 
    Force { Vector3d::Zero() }, 
    VerletPosition { Vector3d::Zero() }, 
    VerletList { } { }

    
MDParticle::MDParticle(Vector3d aPosition, Vector3d aVelocity, double aMass) : 
    MPCParticle(aPosition, aVelocity, aMass), 
    Force { Vector3d::Zero() }, 
    VerletPosition { Vector3d::Zero() }, 
    VerletList { } { }
    