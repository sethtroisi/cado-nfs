/* This must be compiled with gcc -x assembler.
 *
 * It's an outright silly test, in order to make a guess as to whether
 * linalg/bwc/matmul-sliced-asm.S should be compiled or not.  The code does
 * not mean anything, especially not for main ! But cmake insists on having
 * test cases have a main() function...
 */
	.text
.Ltext0:
	.p2align 4,,15
.globl main
	.type	main, @function
main:
        pushq %rbp
        pushq %rbx
        pushq %r15
        pushq %r14
        pushq %r13
        pushq %r12
        movq    %rsi, %rdi
	movzwl	2(%rcx), %eax
	movq	%rdx, %rsi
	movzwl	(%rcx), %edx
	addq	$4, %rcx
	salq	$16, %rax
	movq	%rax, %r9
        leaq    (%rcx,%r9,4),%rbx
        salq    $3, %r8
	movq	%rbx, %rax
        popq %r12
        popq %r13
        popq %r14
        popq %r15
        popq %rbx
        popq %rbp
	ret
	.size	main, .-main
