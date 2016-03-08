#ifndef LAS_BASE_HPP_
#define LAS_BASE_HPP_

/* C++ base classes that are not application-specific. */


/* Base class with private copy-constructor and assignment operator.
   Classes which are not copy-constructible can inherit this with:
   private NonCopyable */
class NonCopyable {
 protected:
   NonCopyable() {}
   ~NonCopyable() {}
 private:
   NonCopyable(const NonCopyable&);
   NonCopyable& operator=(const NonCopyable&);
};

#endif
