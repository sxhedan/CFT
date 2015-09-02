##FT library checks

../lib/x86_64/libFronTier.a: ../lib/x86_64/libfront.a ../lib/x86_64/libutil.a ../lib/x86_64/libintfc.a
	@cd ../lib/x86_64/; \
	rm -f *.o; \
	ar x libintfc.a; \
	ar x libutil.a;\
	ar x libfront.a;\
	ar cru libFronTier.a *.o; \
	ranlib libFronTier.a; \
	rm -rf *.o

FT_libs:
	@cd ../src/front && $(MAKE) ../../lib/x86_64/libfront.a
	@cd ../src/util && $(MAKE) ../../lib/x86_64/libutil.a
	@cd ../src/intfc && $(MAKE) ../../lib/x86_64/libintfc.a
	@cd .. && $(MAKE) incs
	@$(MAKE) ../lib/x86_64/libFronTier.a
