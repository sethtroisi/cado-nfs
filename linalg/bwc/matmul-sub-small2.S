	.text
.Ltext0:
	.p2align 4,,15
.globl matmul_sub_small2
	.type	matmul_sub_small2, @function

/* It's possible to do some stuff with SSE instructions ; however it
 * doesn't pay off, or I'm doing it wrong.
 */
#define xxxSSE

#ifdef  SSE
#define ALIGN_BYTES    16
#define ALIGN_MASK     15
#define NEED_ALLOCA
#else
#define ALIGN_BYTES    8
#define ALIGN_MASK     7
#define xxxNEED_ALLOCA
#endif

matmul_sub_small2:
        pushq %rbp
        pushq %rbx
        pushq %r15
        pushq %r14
        pushq %r13
        pushq %r12
        // dst pointer --> store in rdi
        movq    %rsi, %rdi
        // src pointer --> store in rsi
	movq	%rdx, %rsi
        // rcx is the data pointer.

#ifdef  NEED_ALLOCA
        // we have to do some ugly tricks to allocate an aligned zone
        movq %rsp, %rax
        subq $8, %rsp
        movq %rsp, %rbx
        andq $ALIGN_MASK, %rbx
        subq %rbx, %rsp
        movq %rax, (%rsp)
        subq $ALIGN_BYTES, %rsp
        subq $ALIGN_BYTES, %rsp
#endif

        // read ncoeffs_slice into r9
	movzwq	2(%rcx), %rax
	movzwq	(%rcx), %rdx
	addq	$4, %rcx
	salq	$16, %rax
	movq	%rax, %r9
	orq	%rdx, %r9
        leaq    (%rcx,%r9,2),%rbx
        // rbx : end data pointer.

        // r9 : the number of coeffs.
	je	.L2

        // we now want to align our reads.
        movq    %rcx, %rax
        decq    %rax
        orq     $ALIGN_MASK, %rax
        incq    %rax
        // rax is the first aligned address above rcx.

        // now rbx is not necessarily aligned either. If it is, we keep
        // it as is. If it isn't, then we use the low bits of rbx to
        // keep the information on the ending alignment.
        movq    %rbx,%rdx
        andq    $ALIGN_MASK,%rdx
        je      .L31
        subq    $ALIGN_BYTES,%rbx
.L31:
        cmpq %rcx, %rax
        je .L30
        movzwq  (%rcx), %r10
        addq    $2, %rcx
        movq    %r10, %r8
        andq    $4095, %r8
        shrq    $12, %r10
        leaq    (%rsi,%r10,8), %rsi
        movq    (%rsi), %rdx
	xorq	%rdx, (%rdi,%r8,8)
        jmp .L31
.L30:
        // now rcx is aligned. Good.
.L4:
#ifndef SSE
#if 0
        /* considerably slower than the other version below */
	movq	(%rcx), %rax 
        movq    %rax, %r8
        movq    %rax, %r9
        andq    $0xf000, %r8
        andq    $0x0fff, %r9
        shrq    $16, %rax
        movq    %rax, %r10
        movq    %rax, %r11
        andq    $0xf000, %r10
        andq    $0x0fff, %r11
        shrq    $16, %rax
        movq    %rax, %r12
        movq    %rax, %r13
        andq    $0xf000, %r12
        andq    $0x0fff, %r13
        shrq    $16, %rax
        movq    %rax, %r14
        movq    %rax, %r15
        andq    $0xf000, %r14
        andq    $0x0fff, %r15
        shrq    $9, %r8
        shlq    $3, %r9
        shrq    $9, %r10
        shlq    $3, %r11
        shrq    $9, %r12
        shlq    $3, %r13
        shrq    $9, %r14
        shlq    $3, %r15
#else
        movzwq  (%rcx),%r8
        movq    %r8, %r9
        andq    $0xf000, %r8
        andq    $0x0fff, %r9
        shrq    $9, %r8
        shlq    $3, %r9
        movzwq  2(%rcx),%r10
        movq    %r10, %r11
        andq    $0xf000, %r10
        andq    $0x0fff, %r11
        shrq    $9, %r10
        shlq    $3, %r11
        movzwq  4(%rcx),%r12
        movq    %r12, %r13
        andq    $0xf000, %r12
        andq    $0x0fff, %r13
        shrq    $9, %r12
        shlq    $3, %r13
        movzwq  6(%rcx),%r14
        movq    %r14, %r15
        andq    $0xf000, %r14
        andq    $0x0fff, %r15
        shrq    $9, %r14
        shlq    $3, %r15
#endif


	movq	(%rsi,%r8), %rdx
	xorq	%rdx,(%rdi,%r9)
        addq    %r10, %r8

	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r11)
        addq    %r12, %r8

	movq	(%rsi,%r8), %rdx
	xorq	%rdx, (%rdi,%r13)
        addq    %r14, %r8

	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r15)
        addq    %r8, %rsi

	addq	$8, %rcx
#elif ALIGN_BYTES == 16
        movdqa  (%rcx), %xmm0
        movdqa  %xmm0, %xmm1
        psrlw   $12, %xmm0
        psllw   $3, %xmm0
        psllw   $4, %xmm1
        psrlw   $1, %xmm1
        movdqa  %xmm0, (%rsp)
        movdqa  %xmm1, 16(%rsp)

        movzwq  0(%rsp),%r8
        movzwq  16(%rsp),%r9
        movzwq  2(%rsp),%r10
        movzwq  18(%rsp),%r11
        movzwq  4(%rsp),%r12
        movzwq  20(%rsp),%r13
        movzwq  6(%rsp),%r14
        movzwq  22(%rsp),%r15

	movq	(%rsi,%r8), %rdx
	xorq	%rdx,(%rdi,%r9)
        addq    %r10, %r8

	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r11)
        addq    %r12, %r8

	movq	(%rsi,%r8), %rdx
	xorq	%rdx, (%rdi,%r13)
        addq    %r14, %r8

	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r15)
        addq    %r8, %rsi

        movzwq  8(%rsp),%r8
        movzwq  24(%rsp),%r9
        movzwq  10(%rsp),%r10
        movzwq  26(%rsp),%r11
        movzwq  12(%rsp),%r12
        movzwq  28(%rsp),%r13
        movzwq  14(%rsp),%r14
        movzwq  30(%rsp),%r15

	movq	(%rsi,%r8), %rdx
	xorq	%rdx,(%rdi,%r9)
        addq    %r10, %r8

	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r11)
        addq    %r12, %r8

	movq	(%rsi,%r8), %rdx
	xorq	%rdx, (%rdi,%r13)
        addq    %r14, %r8

	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r15)
        addq    %r8, %rsi

	addq	$16, %rcx
#endif

	cmpq	%rcx, %rbx
	ja	.L4

        je      .L40
        addq    $ALIGN_BYTES,%rbx
.L41:
        cmpq %rcx, %rbx
        je .L40
        movzwq  (%rcx), %r10
        addq    $2, %rcx
        movq    %r10, %r8
        andq    $4095, %r8
        shrq    $12, %r10
        leaq    (%rsi,%r10,8), %rsi
        movq    (%rsi), %rdx
	xorq	%rdx, (%rdi,%r8,8)
        jmp .L41
.L40:
.L2:
        // leave.
#ifdef  NEED_ALLOCA
        addq $ALIGN_BYTES, %rsp
        addq $ALIGN_BYTES, %rsp
        movq (%rsp), %rax
        movq %rax, %rsp
#endif

	movq	%rcx, %rax
        popq %r12
        popq %r13
        popq %r14
        popq %r15
        popq %rbx
        popq %rbp
	ret
.LFE52:
	.size	matmul_sub_small2, .-matmul_sub_small2
