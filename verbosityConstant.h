#pragma once

//we're adding a global state, so proceed with caution.
//e.g. Update()s in the verbosity might update the itk pipeline and dodge an otherwise present error

namespace VerbosityConstant {
  enum Verbosity {NONE, LOW, MEDIUM, HIGH};
  static Verbosity verbosity;
}

//using singleton pattern is an overkill. We allow overwriting but well.

//class VerbosityConstant
//{
//public:

//    enum Verbosity {NONE, LOW, MEDIUM, HIGH};

//    // Some accessor functions for the class, itself
//    Verbosity GetValue() const
//    {return m_verbosity;}

//    // We could hide this if we want and allow access only to one friend class/method
//    void SetValue(const Verbosity &verbosity)
//    {m_verbosity = verbosity;}

//    // The magic function, which allows access to the class from anywhere
//    // To get the value of the instance of the class, call:
//    //     VerbosityConstant::Instance().GetValue();
//    static VerbosityConstant &Instance()
//    {
//        // This line only runs once, thus creating the only instance in existence
//        static VerbosityConstant *instance = new VerbosityConstant;
//        // dereferencing the variable here, saves the caller from having to use
//        // the arrow operator, and removes temptation to try and delete the
//        // returned instance.
//        return *instance; // always returns the same instance
//    }

//private:
//    // We need to make some given functions private to finish the definition of the singleton
//    VerbosityConstant(){} // default constructor available only to members or friends of this class

//    // Note that the next two functions are not given bodies, thus any attempt
//    // to call them implicitly will return as compiler errors. This prevents
//    // accidental copying of the only instance of the class.
//    VerbosityConstant(const VerbosityConstant &old); // disallow copy constructor
//    const VerbosityConstant &operator=(const VerbosityConstant &old); //disallow assignment operator

//    // Note that although this should be allowed,
//    // some compilers may not implement private destructors
//    // This prevents others from deleting our one single instance, which was otherwise created on the heap
//    ~VerbosityConstant(){}
//private: // private data for an instance of this class
//    Verbosity m_verbosity;
//};

//class VerbosityConstant
//{
//public:

//    enum Verbosity {NONE, LOW, MEDIUM, HIGH};

//    static int get()
//    {
//        return someConstant;
//    }

//private:
//    friend int main(); // this should probably not be `main` if we want to plug it in Espina

//    static void set(const Verbosity& value)
//    {
//        someConstant = value;
//    }

//    VerbosityConstant() : someConstant{LOW}{
//    } // default constructor available only to members or friends of this class

//    VerbosityConstant(const VerbosityConstant &old); // disallow copy constructor
//    //const VerbosityConstant &operator=(const VerbosityConstant &old); //disallow assignment operator

//    // This prevents others from deleting our one single instance, which was otherwise created on the heap
//    ~VerbosityConstant(){}

//  private:
//    static Verbosity const someConstant;
//};

